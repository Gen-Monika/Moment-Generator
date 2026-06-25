#pragma once

#include "MomentEngine.h"

#include <QMainWindow>
#include <QSet>

#include <optional>

class QCheckBox;
class QLabel;
class QPlainTextEdit;
class QPushButton;
class QSpinBox;
class QTableWidget;
class QWebEngineView;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget* parent = nullptr);

private:
    void buildUi();
    void connectSignals();
    void regenerate();
    void copyLatexToClipboard();
    void populateTable(const QVector<MomentTerm>& terms);
    QString renderHtml(const QVector<FormulaSection>& sections) const;

    QSpinBox* m_orderSpin = nullptr;
    QSpinBox* m_sampleSizeSpin = nullptr;
    QCheckBox* m_symbolicNCheck = nullptr;
    QPushButton* m_generateButton = nullptr;
    QPushButton* m_copyButton = nullptr;
    QTableWidget* m_partitionTable = nullptr;
    QWebEngineView* m_preview = nullptr;
    QPlainTextEdit* m_latexEdit = nullptr;
    QLabel* m_statusLabel = nullptr;
    QSet<int> m_warnedOnDemandOrders;
};


