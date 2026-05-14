import sys
import os
import io
import base64
import pickle
import threading
import webbrowser
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from flask import Flask, send_from_directory, request, jsonify

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import project as ml

frontend_dir = os.path.abspath(os.path.dirname(__file__))
app = Flask(__name__, static_folder=frontend_dir, static_url_path="")


def _problem_to_json(prob):
    blob = base64.b64encode(pickle.dumps(prob)).decode()
    return {
        "kind": prob.kind,
        "prompt": prob.prompt,
        "prompt_latex": getattr(prob, "prompt_latex", prob.prompt),
        "_blob": blob,
    }


def _json_to_problem(data):
    return pickle.loads(base64.b64decode(data["_blob"]))


def _make_plot_data(report):
    try:
        obj = report.answer_expr
        if obj is None:
            return None
        x = sp.Symbol("x")
        free = obj.free_symbols if hasattr(obj, "free_symbols") else set()
        if free and free != {x}:
            return None
        expr = obj.rhs if isinstance(obj, sp.Eq) else obj
        f = sp.lambdify(x, expr, modules=["numpy"])
        xs = np.linspace(-5, 5, 400)
        ys = f(xs).astype(float)
        ys[np.abs(ys) > 100] = np.nan
        return {"x": xs.tolist(), "y": [None if np.isnan(v) else v for v in ys]}
    except Exception:
        return None


def _render_plot_png(x_data, y_data, kind):
    fig, ax = plt.subplots(figsize=(6, 4))
    xs = np.array(x_data)
    ys = np.array([np.nan if v is None else v for v in y_data])
    ax.plot(xs, ys, linewidth=2)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(kind.replace("_", " "))
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5)
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=100, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()


def _report_to_json(report, want_plot_data):
    steps = []
    for s in report.steps:
        steps.append({
            "title": s.title,
            "math_latex": s.math,
            "explanation": s.explanation,
        })

    try:
        answer_latex = sp.latex(report.answer_expr)
    except Exception:
        answer_latex = report.final_answer

    plot_data = None
    if want_plot_data:
        plot_data = _make_plot_data(report)

    return {
        "kind": report.kind,
        "final_answer_latex": answer_latex,
        "steps": steps,
        "verified": report.verified,
        "verify_msg": report.verify_msg,
        "plot_data": plot_data,
        "warnings": report.warnings,
    }


@app.route("/")
def index():
    return send_from_directory(frontend_dir, "index.html")


@app.route("/api/kinds")
def api_kinds():
    return jsonify({"kinds": ml.list_kinds()})


@app.route("/api/generate", methods=["POST"])
def api_generate():
    try:
        body = request.get_json()
        kind = body.get("kind", "ode1_linear")
        difficulty = body.get("difficulty", "easy")
        with_ics = body.get("with_ics", True)
        prob = ml.generate(kind, difficulty=difficulty, with_ics=with_ics)
        return jsonify({"problem": _problem_to_json(prob)})
    except ml.TutorError as e:
        return jsonify({"error": str(e)}), 400


@app.route("/api/solve", methods=["POST"])
def api_solve():
    try:
        body = request.get_json()
        prob = _json_to_problem(body["problem"])
        want_plot_data = body.get("want_plot_data", True)
        report = ml.solve(prob)
        result = _report_to_json(report, want_plot_data)
        return jsonify({"report": result})
    except ml.TutorError as e:
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/plot", methods=["POST"])
def api_plot():
    try:
        body = request.get_json()
        plot_data = body["plot_data"]
        kind = body.get("kind", "")
        png_b64 = _render_plot_png(plot_data["x"], plot_data["y"], kind)
        return jsonify({"png_base64": png_b64})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


if __name__ == "__main__":
    threading.Timer(1.0, lambda: webbrowser.open("http://127.0.0.1:8080")).start()
    app.run(debug=True, use_reloader=False)
