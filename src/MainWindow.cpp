#include "MainWindow.h"

#include <QApplication>
#include <QCheckBox>
#include <QClipboard>
#include <QHeaderView>
#include <QLabel>
#include <QPlainTextEdit>
#include <QPushButton>
#include <QSpinBox>
#include <QSplitter>
#include <QTableWidget>
#include <QVBoxLayout>
#include <QWebEngineView>

MainWindow::MainWindow(QWidget* parent)
    : QMainWindow(parent)
{
    buildUi();
    connectSignals();
    regenerate();
}

void MainWindow::buildUi()
{
    setWindowTitle(QStringLiteral("Moment Generator"));
    resize(1280, 820);

    auto* central = new QWidget(this);
    auto* rootLayout = new QVBoxLayout(central);
    rootLayout->setContentsMargins(12, 12, 12, 12);

    auto* splitter = new QSplitter(Qt::Horizontal, central);

    auto* leftPanel = new QWidget(splitter);
    auto* leftLayout = new QVBoxLayout(leftPanel);
    leftLayout->setContentsMargins(0, 0, 8, 0);

    auto* controls = new QWidget(leftPanel);
    auto* controlsLayout = new QVBoxLayout(controls);
    controlsLayout->setContentsMargins(0, 0, 0, 0);

    auto* orderLabel = new QLabel(QStringLiteral("Moment order k"), controls);
    m_orderSpin = new QSpinBox(controls);
    m_orderSpin->setRange(2, 20);
    m_orderSpin->setValue(4);

    auto* sampleLabel = new QLabel(QStringLiteral("Sample size n"), controls);
    m_sampleSizeSpin = new QSpinBox(controls);
    m_sampleSizeSpin->setRange(2, 1000000);
    m_sampleSizeSpin->setValue(10);

    m_symbolicNCheck = new QCheckBox(QStringLiteral("Keep n symbolic"), controls);
    m_symbolicNCheck->setChecked(true);

    m_generateButton = new QPushButton(QStringLiteral("Generate"), controls);
    m_copyButton = new QPushButton(QStringLiteral("Copy LaTeX"), controls);

    controlsLayout->addWidget(orderLabel);
    controlsLayout->addWidget(m_orderSpin);
    controlsLayout->addWidget(sampleLabel);
    controlsLayout->addWidget(m_sampleSizeSpin);
    controlsLayout->addWidget(m_symbolicNCheck);
    controlsLayout->addWidget(m_generateButton);
    controlsLayout->addWidget(m_copyButton);

    m_partitionTable = new QTableWidget(leftPanel);
    m_partitionTable->setColumnCount(4);
    m_partitionTable->setHorizontalHeaderLabels({
        QStringLiteral("Partition"),
        QStringLiteral("Product"),
        QStringLiteral("Coefficient"),
        QStringLiteral("Numeric")
    });
    m_partitionTable->horizontalHeader()->setStretchLastSection(true);
    m_partitionTable->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    m_partitionTable->verticalHeader()->setVisible(false);
    m_partitionTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_partitionTable->setSelectionBehavior(QAbstractItemView::SelectRows);

    leftLayout->addWidget(controls);
    leftLayout->addWidget(new QLabel(QStringLiteral("Integer partitions without 1"), leftPanel));
    leftLayout->addWidget(m_partitionTable, 1);

    auto* rightPanel = new QWidget(splitter);
    auto* rightLayout = new QVBoxLayout(rightPanel);
    rightLayout->setContentsMargins(8, 0, 0, 0);

    m_preview = new QWebEngineView(rightPanel);
    m_latexEdit = new QPlainTextEdit(rightPanel);
    m_latexEdit->setPlaceholderText(QStringLiteral("Generated LaTeX will appear here."));
    m_latexEdit->setMinimumHeight(180);

    rightLayout->addWidget(new QLabel(QStringLiteral("KaTeX preview"), rightPanel));
    rightLayout->addWidget(m_preview, 2);
    rightLayout->addWidget(new QLabel(QStringLiteral("LaTeX source"), rightPanel));
    rightLayout->addWidget(m_latexEdit, 1);

    splitter->addWidget(leftPanel);
    splitter->addWidget(rightPanel);
    splitter->setStretchFactor(0, 0);
    splitter->setStretchFactor(1, 1);

    m_statusLabel = new QLabel(central);

    rootLayout->addWidget(splitter, 1);
    rootLayout->addWidget(m_statusLabel);
    setCentralWidget(central);
}

void MainWindow::connectSignals()
{
    connect(m_generateButton, &QPushButton::clicked, this, &MainWindow::regenerate);
    connect(m_copyButton, &QPushButton::clicked, this, &MainWindow::copyLatexToClipboard);
    connect(m_orderSpin, qOverload<int>(&QSpinBox::valueChanged), this, &MainWindow::regenerate);
    connect(m_sampleSizeSpin, qOverload<int>(&QSpinBox::valueChanged), this, &MainWindow::regenerate);
    connect(m_symbolicNCheck, &QCheckBox::toggled, this, &MainWindow::regenerate);
}

void MainWindow::regenerate()
{
    const int order = m_orderSpin->value();
    const std::optional<int> sampleSize = m_symbolicNCheck->isChecked()
        ? std::nullopt
        : std::optional<int>(m_sampleSizeSpin->value());

    const QVector<MomentTerm> terms = MomentEngine::terms(order, sampleSize);
    const QString latex = MomentEngine::fullLatex(order, terms, sampleSize);

    populateTable(terms);
    m_latexEdit->setPlainText(latex);
    m_preview->setHtml(renderHtml(latex), QUrl(QStringLiteral("https://cdn.jsdelivr.net/")));
    m_statusLabel->setText(QStringLiteral("Generated %1 partition term(s) for k=%2.").arg(terms.size()).arg(order));
}

void MainWindow::copyLatexToClipboard()
{
    QApplication::clipboard()->setText(m_latexEdit->toPlainText());
    m_statusLabel->setText(QStringLiteral("LaTeX copied to clipboard."));
}

void MainWindow::populateTable(const QVector<MomentTerm>& terms)
{
    m_partitionTable->setRowCount(terms.size());
    for (int row = 0; row < terms.size(); ++row) {
        const MomentTerm& term = terms[row];
        m_partitionTable->setItem(row, 0, new QTableWidgetItem(term.partition.labelLatex()));
        m_partitionTable->setItem(row, 1, new QTableWidgetItem(term.productLatex));
        m_partitionTable->setItem(row, 2, new QTableWidgetItem(term.coefficientLatex));
        m_partitionTable->setItem(row, 3, new QTableWidgetItem(
            term.numericCoefficient.has_value() ? QString::number(static_cast<double>(*term.numericCoefficient), 'g', 10)
                                                : QStringLiteral("-")));
    }
}

QString MainWindow::renderHtml(const QString& latex) const
{
    const QString escapedLatex = latex.toHtmlEscaped();
    return QString::fromUtf8(R"HTML(<!doctype html>
<html>
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.css">
  <script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/katex.min.js"></script>
  <script defer src="https://cdn.jsdelivr.net/npm/katex@0.16.11/dist/contrib/auto-render.min.js"></script>
  <style>
    body {
      margin: 0;
      padding: 24px;
      font-family: "Segoe UI", "Microsoft YaHei", sans-serif;
      color: #202124;
      background: #ffffff;
    }
    .formula {
      overflow-x: auto;
      padding: 18px;
      border: 1px solid #d8dee4;
      background: #f8f8f8;
    }
    .hint {
      margin-top: 14px;
      color: #666;
      font-size: 13px;
    }
  </style>
</head>
<body>
  <div class="formula">\[%1\]</div>
  <div class="hint">KaTeX is loaded from jsDelivr in this first scaffold. Offline assets can be vendored later.</div>
  <script>
    document.addEventListener("DOMContentLoaded", function() {
      renderMathInElement(document.body, {
        delimiters: [
          {left: "\\[", right: "\\]", display: true},
          {left: "$$", right: "$$", display: true},
          {left: "\\(", right: "\\)", display: false}
        ],
        throwOnError: false
      });
    });
  </script>
</body>
</html>)HTML").arg(escapedLatex);
}


