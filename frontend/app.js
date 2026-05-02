// I am so so bad at all frontend languages so I used help from my friend

let currentProblem = null;

async function loadKinds() {
  const res = await fetch("/api/kinds");
  const data = await res.json();
  const select = document.getElementById("kindSelect");
  for (const kind of data.kinds) {
    const opt = document.createElement("option");
    opt.value = kind;
    opt.textContent = kind.replace(/_/g, " ");
    select.appendChild(opt);
  }
}

async function generate() {
  const kind = document.getElementById("kindSelect").value;
  const difficulty = document.getElementById("difficultySelect").value;
  const withIcs = document.getElementById("icsToggle").checked;

  const body = { kind, difficulty, with_ics: withIcs };

  document.getElementById("generateBtn").disabled = true;
  const res = await fetch("/api/generate", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  document.getElementById("generateBtn").disabled = false;
  const data = await res.json();
  if (data.error) { alert(data.error); return; }

  currentProblem = data.problem;

  var card = document.getElementById("problemPrompt");
  card.className = "card";
  card.innerHTML = "<h2>Problem</h2><p>$$" + currentProblem.prompt_latex + "$$</p>";
  MathJax.typesetPromise([card]);

  document.getElementById("solveBtn").disabled = false;

  document.getElementById("stepsContainer").className = "hidden";
  document.getElementById("answerCard").className = "card placeholder";
  document.getElementById("answerCard").innerHTML = '<h2>Result</h2><p class="placeholder-text">Solve a problem to see the answer</p>';
  document.getElementById("plotCard").className = "card hidden";
}

async function solve() {
  if (!currentProblem) return;

  document.getElementById("spinner").classList.remove("hidden");
  document.getElementById("solveBtn").disabled = true;

  const res = await fetch("/api/solve", {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      problem: currentProblem,
      want_plot_data: true,
    }),
  });
  const data = await res.json();

  document.getElementById("spinner").classList.add("hidden");
  document.getElementById("solveBtn").disabled = false;

  if (data.error) { alert(data.error); return; }

  var report = data.report;

  var container = document.getElementById("stepsContainer");
  var html = "<h2>Solution Steps</h2>";
  report.steps.forEach(function(step, i) {
    html += '<div class="step-card" onclick="this.querySelector(\'.step-body\').classList.toggle(\'collapsed\')">'
      + "<strong>Step " + (i + 1) + ": " + step.title + "</strong>"
      + '<div class="step-body">'
      + "<p>$$" + step.math_latex + "$$</p>"
      + '<p class="explanation">' + step.explanation + "</p>"
      + "</div></div>";
  });
  container.innerHTML = html;
  container.classList.remove("hidden");

  var answerCard = document.getElementById("answerCard");
  var badge = report.verified
    ? '<span id="verifiedBadge" class="badge pass">Verified</span>'
    : '<span id="verifiedBadge" class="badge fail">Unverified</span>';
  var answerHtml = "<h2>Final Answer " + badge + "</h2>"
    + '<div id="finalAnswer" class="answer-expr">$$' + report.final_answer_latex + "$$</div>";
  if (report.warnings && report.warnings.length > 0) {
    answerHtml += '<p class="warnings">Warnings: ' + report.warnings.join(", ") + "</p>";
  }
  answerCard.className = "card";
  answerCard.innerHTML = answerHtml;

  if (report.plot_data) {
    var plotRes = await fetch("/api/plot", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ plot_data: report.plot_data, kind: report.kind }),
    });
    var plotJson = await plotRes.json();
    if (plotJson.png_base64) {
      document.getElementById("plotImg").src = "data:image/png;base64," + plotJson.png_base64;
      document.getElementById("plotCard").className = "card";
    }
  } else {
    document.getElementById("plotCard").className = "card hidden";
  }

  MathJax.typesetPromise([container, answerCard]);
}

document.getElementById("generateBtn").addEventListener("click", generate);
document.getElementById("solveBtn").addEventListener("click", solve);

document.getElementById("plotImg").addEventListener("click", function() {
  this.classList.toggle("zoomed");
});

loadKinds();
