#include "MomentEngine.h"

#include <QStringList>
#include <QtMath>

#include <algorithm>
#include <cmath>

QMap<int, int> Partition::multiplicities() const
{
    QMap<int, int> counts;
    for (int part : parts) {
        counts[part] += 1;
    }
    return counts;
}

int Partition::blockCount() const
{
    return parts.size();
}

QString Partition::labelLatex() const
{
    QStringList labels;
    labels.reserve(parts.size());
    for (int part : parts) {
        labels << QString::number(part);
    }
    return QStringLiteral("(%1)").arg(labels.join(QStringLiteral(",")));
}

QString Partition::productLatex() const
{
    const QMap<int, int> counts = multiplicities();
    QStringList factors;
    for (auto it = counts.cbegin(); it != counts.cend(); ++it) {
        if (it.value() == 1) {
            factors << QStringLiteral("v_{%1}").arg(it.key());
        } else {
            factors << QStringLiteral("v_{%1}^{%2}").arg(it.key()).arg(it.value());
        }
    }
    std::reverse(factors.begin(), factors.end());
    return factors.join(QStringLiteral(" "));
}

QVector<Partition> MomentEngine::partitionsWithoutOnes(int order)
{
    QVector<Partition> out;
    QVector<int> current;
    generatePartitions(order, order, current, out);
    return out;
}

QVector<MomentTerm> MomentEngine::terms(int order, std::optional<int> sampleSize)
{
    QVector<MomentTerm> result;
    const QVector<Partition> partitions = partitionsWithoutOnes(order);
    result.reserve(partitions.size());

    for (const Partition& partition : partitions) {
        MomentTerm term;
        term.partition = partition;
        term.coefficientSymbolLatex = QStringLiteral("c_{%1}(n)").arg(partition.labelLatex());
        term.coefficientLatex = coefficientLatex(order, partition);
        term.productLatex = partition.productLatex();
        term.symbolicTermLatex = QStringLiteral("%1\\, %2").arg(term.coefficientSymbolLatex, term.productLatex);

        if (sampleSize.has_value()) {
            term.numericCoefficient = numericCoefficient(order, partition, *sampleSize);
            term.numericTermLatex = QStringLiteral("%1\\,%2").arg(formatNumber(*term.numericCoefficient), term.productLatex);
        }

        result << term;
    }

    return result;
}

QVector<FormulaSection> MomentEngine::formulaSections(int order, const QVector<MomentTerm>& terms, std::optional<int> sampleSize)
{
    QVector<FormulaSection> sections;
    sections << FormulaSection{QStringLiteral("Definitions"), definitionsLatex(order)};
    sections << FormulaSection{QStringLiteral("Forward relation"), forwardLatex(order, terms, sampleSize)};
    sections << FormulaSection{QStringLiteral("Inverse estimator"), inverseLatex(order, terms)};

    const QString applications = applicationsLatex(order);
    if (!applications.isEmpty()) {
        sections << FormulaSection{QStringLiteral("Applications"), applications};
    }

    return sections;
}

QString MomentEngine::fullLatex(int order, const QVector<MomentTerm>& terms, std::optional<int> sampleSize)
{
    const QVector<FormulaSection> sections = formulaSections(order, terms, sampleSize);
    QStringList blocks;
    for (const FormulaSection& section : sections) {
        blocks << QStringLiteral("% %1\n\\[\n%2\n\\]")
                      .arg(section.title, section.latex);
    }
    return blocks.join(QStringLiteral("\n\n"));
}

void MomentEngine::generatePartitions(int remaining, int maxPart, QVector<int>& current, QVector<Partition>& out)
{
    if (remaining == 0) {
        out << Partition{current};
        return;
    }

    for (int part = std::min(maxPart, remaining); part >= 2; --part) {
        current << part;
        generatePartitions(remaining - part, part, current, out);
        current.removeLast();
    }
}

qint64 MomentEngine::factorial(int value)
{
    qint64 result = 1;
    for (int i = 2; i <= value; ++i) {
        result *= i;
    }
    return result;
}

qint64 MomentEngine::multinomialBlockFactor(int order, const QMap<int, int>& multiplicities)
{
    qint64 denominator = 1;
    for (auto it = multiplicities.cbegin(); it != multiplicities.cend(); ++it) {
        for (int count = 0; count < it.value(); ++count) {
            denominator *= factorial(it.key());
        }
        denominator *= factorial(it.value());
    }
    return factorial(order) / denominator;
}

QString MomentEngine::coefficientLatex(int order, const Partition& partition)
{
    const QMap<int, int> counts = partition.multiplicities();
    const int blocks = partition.blockCount();
    const qint64 factor = multinomialBlockFactor(order, counts);
    const QString falling = fallingFactorLatex(1, blocks - 1);
    const QString bracket = bracketLatex(order, counts, blocks);

    QString numerator;
    if (factor != 1) {
        numerator += QString::number(factor);
    }
    if (falling != QStringLiteral("1")) {
        numerator += falling;
    }
    numerator += QStringLiteral("\\left[%1\\right]").arg(bracket);

    return QStringLiteral("\\frac{%1}{n^{%2}}").arg(numerator).arg(order);
}

QString MomentEngine::bracketLatex(int order, const QMap<int, int>& multiplicities, int blockCount)
{
    QStringList pieces;

    for (auto it = multiplicities.cbegin(); it != multiplicities.cend(); ++it) {
        const int q = it.key();
        const int count = it.value();
        const bool positive = ((order - q) % 2 == 0);

        QString body;
        if (count != 1) {
            body += QString::number(count);
        }
        body += QStringLiteral("(n-1)");
        if (q != 1) {
            body += QStringLiteral("^{%1}").arg(q);
        }

        pieces << (positive ? QStringLiteral("+") : QStringLiteral("-")) + body;
    }

    const bool tailPositive = (order % 2 == 0);
    pieces << (tailPositive ? QStringLiteral("+") : QStringLiteral("-"))
              + QStringLiteral("(n-%1)").arg(blockCount);

    QString result;
    for (const QString& piece : pieces) {
        if (result.isEmpty()) {
            result = piece.startsWith(QLatin1Char('+')) ? piece.mid(1) : piece;
        } else {
            result += piece.startsWith(QLatin1Char('+')) ? QStringLiteral("+") + piece.mid(1) : piece;
        }
    }
    return result;
}

QString MomentEngine::fallingFactorLatex(int start, int count)
{
    if (count <= 0) {
        return QStringLiteral("1");
    }

    QStringList factors;
    for (int i = 0; i < count; ++i) {
        const int offset = start + i;
        factors << QStringLiteral("(n-%1)").arg(offset);
    }
    return factors.join(QString());
}

QString MomentEngine::definitionsLatex(int order)
{
    return QStringLiteral(
        "\\begin{aligned}\n"
        "M_{%1} &= \\frac{1}{n}\\sum_{i=1}^{n}(X_i-\\bar X)^{%1},\\\\[0.35em]\n"
        "v_r &= \\mathbb{E}[(X-\\mu)^r],\\qquad v_1=0,\\\\[0.35em]\n"
        "\\widehat v_r &\\text{ denotes the recursively corrected estimator of } v_r,\\\\[0.35em]\n"
        "\\mathbb{S},\\mathbb{K},\\mathbb{EK} &\\text{ denote skewness, kurtosis, and excess kurtosis.}\n"
        "\\end{aligned}")
        .arg(order);
}

QString MomentEngine::forwardLatex(int order, const QVector<MomentTerm>& terms, std::optional<int> sampleSize)
{
    QStringList mainTerms;
    QStringList coefficientLines;

    for (const MomentTerm& term : terms) {
        mainTerms << term.symbolicTermLatex;
        coefficientLines << QStringLiteral("%1 &= %2").arg(term.coefficientSymbolLatex, term.coefficientLatex);
    }

    QString latex;
    latex += QStringLiteral("\\begin{aligned}\n");
    latex += QStringLiteral("\\mathbb{E}[M_{%1}] &= %2\\\\[0.8em]\n")
                 .arg(order)
                 .arg(joinTerms(mainTerms));
    latex += coefficientLines.join(QStringLiteral("\\\\[0.45em]\n"));

    if (sampleSize.has_value()) {
        QStringList numericTerms;
        for (const MomentTerm& term : terms) {
            numericTerms << term.numericTermLatex;
        }
        latex += QStringLiteral("\\\\[0.8em]\n");
        latex += QStringLiteral("\\mathbb{E}[M_{%1}]\\big|_{n=%2} &= %3")
                     .arg(order)
                     .arg(*sampleSize)
                     .arg(joinTerms(numericTerms));
    }

    latex += QStringLiteral("\n\\end{aligned}");
    return latex;
}

QString MomentEngine::inverseLatex(int order, const QVector<MomentTerm>& terms)
{
    const MomentTerm* lead = leadingTerm(terms, order);
    if (lead == nullptr) {
        return QStringLiteral("\\text{No leading partition was found.}");
    }

    QStringList lowerExpectationTerms;
    QStringList lowerEstimatorTerms;
    for (const MomentTerm& term : terms) {
        if (&term == lead) {
            continue;
        }
        lowerExpectationTerms << term.symbolicTermLatex;
        lowerEstimatorTerms << QStringLiteral("%1\\, %2")
                                   .arg(term.coefficientSymbolLatex, estimatorProductLatex(term.partition));
    }

    const QString expectationNumerator = lowerExpectationTerms.isEmpty()
        ? QStringLiteral("\\mathbb{E}[M_{%1}]").arg(order)
        : QStringLiteral("\\mathbb{E}[M_{%1}] - \\left(%2\\right)")
              .arg(order)
              .arg(joinTerms(lowerExpectationTerms));

    const QString estimatorNumerator = lowerEstimatorTerms.isEmpty()
        ? QStringLiteral("M_{%1}").arg(order)
        : QStringLiteral("M_{%1} - \\left(%2\\right)")
              .arg(order)
              .arg(joinTerms(lowerEstimatorTerms));

    return QStringLiteral(
        "\\begin{aligned}\n"
        "v_{%1} &= \\frac{%2}{%3},\\\\[0.8em]\n"
        "\\widehat v_{%1} &= \\frac{%4}{%3}.\n"
        "\\end{aligned}")
        .arg(order)
        .arg(expectationNumerator)
        .arg(lead->coefficientSymbolLatex)
        .arg(estimatorNumerator);
}

QString MomentEngine::applicationsLatex(int order)
{
    if (order == 3) {
        return QStringLiteral(
            "\\begin{aligned}\n"
            "\\mathbb{S}(X) &= \\frac{v_3}{v_2^{3/2}},\\qquad "
            "\\mathbb{S}_{M}(X)=\\frac{M_3}{M_2^{3/2}},\\\\[0.45em]\n"
            "\\widehat v_2 &= \\frac{n}{n-1}M_2,\\qquad "
            "\\widehat v_3=\\frac{n^2}{(n-1)(n-2)}M_3,\\\\[0.45em]\n"
            "\\widehat{\\mathbb{S}}_{\\mathrm{plug\\text{-}in}}(X) "
            "&= \\frac{\\widehat v_3}{\\widehat v_2^{3/2}}\\\\[0.35em]\n"
            "&= \\frac{\\sqrt{n(n-1)}}{n-2}\\,"
            "\\frac{M_3}{M_2^{3/2}},\\qquad n>2.\n"
            "\\end{aligned}");
    }

    if (order == 4) {
        return QStringLiteral(
            "\\begin{aligned}\n"
            "\\mathbb{K}(X) &= \\frac{v_4}{v_2^2},\\qquad "
            "\\mathbb{EK}(X)=\\mathbb{K}(X)-3=\\frac{v_4-3v_2^2}{v_2^2},\\\\[0.35em]\n"
            "\\mathbb{K}_{M}(X) &= \\frac{M_4}{M_2^2},\\qquad "
            "\\mathbb{EK}_{M}(X)=\\mathbb{K}_{M}(X)-3,\\\\[0.55em]\n"
            "\\widehat v_2 &= \\frac{n}{n-1}M_2,\\\\[0.35em]\n"
            "\\widehat v_4 &= \\frac{n(n^2-2n+3)}{(n-1)(n-2)(n-3)}M_4"
            "-\\frac{3n(2n-3)}{(n-1)(n-2)(n-3)}M_2^2,\\\\[0.55em]\n"
            "\\widehat{\\mathbb{K}}_{\\mathrm{plug\\text{-}in}}(X) "
            "&= \\frac{\\widehat v_4}{\\widehat v_2^2}\\\\[0.35em]\n"
            "&= \\frac{(n-1)(n^2-2n+3)}{n(n-2)(n-3)}\\,\\frac{M_4}{M_2^2}"
            "-\\frac{3(n-1)(2n-3)}{n(n-2)(n-3)},\\\\[0.55em]\n"
            "\\widehat{\\mathbb{EK}}_{\\mathrm{plug\\text{-}in}}(X) "
            "&= \\widehat{\\mathbb{K}}_{\\mathrm{plug\\text{-}in}}(X)-3\\\\[0.35em]\n"
            "&= \\frac{(n-1)(n^2-2n+3)}{n(n-2)(n-3)}\\,\\frac{M_4}{M_2^2}"
            "-\\frac{3(n^3-3n^2+n+3)}{n(n-2)(n-3)},\\qquad n>3,\\\\[0.55em]\n"
            "\\widehat\\kappa_4 &= \\frac{n^2\\left[(n+1)M_4-3(n-1)M_2^2\\right]}{(n-1)(n-2)(n-3)},\\\\[0.35em]\n"
            "\\widehat{\\mathbb{EK}}_{\\kappa,\\mathrm{plug\\text{-}in}}(X) "
            "&= \\frac{\\widehat\\kappa_4}{\\widehat v_2^2}"
            "=\\frac{n-1}{(n-2)(n-3)}\\left[(n+1)\\frac{M_4}{M_2^2}-3(n-1)\\right].\n"
            "\\end{aligned}");
    }

    return QString();
}

QString MomentEngine::estimatorProductLatex(const Partition& partition)
{
    const QMap<int, int> counts = partition.multiplicities();
    QStringList factors;
    for (auto it = counts.cbegin(); it != counts.cend(); ++it) {
        if (it.value() == 1) {
            factors << QStringLiteral("\\widehat v_{%1}").arg(it.key());
        } else {
            factors << QStringLiteral("\\widehat v_{%1}^{%2}").arg(it.key()).arg(it.value());
        }
    }
    std::reverse(factors.begin(), factors.end());
    return factors.join(QStringLiteral(" "));
}

QString MomentEngine::joinTerms(const QStringList& terms)
{
    QString result;
    for (QString term : terms) {
        term = term.trimmed();
        if (term.isEmpty()) {
            continue;
        }

        if (result.isEmpty()) {
            result = term;
        } else if (term.startsWith(QLatin1Char('-'))) {
            result += QStringLiteral(" - ") + term.mid(1).trimmed();
        } else {
            result += QStringLiteral(" + ") + term;
        }
    }
    return result;
}

const MomentTerm* MomentEngine::leadingTerm(const QVector<MomentTerm>& terms, int order)
{
    for (const MomentTerm& term : terms) {
        if (term.partition.parts.size() == 1 && term.partition.parts.first() == order) {
            return &term;
        }
    }
    return nullptr;
}

long double MomentEngine::numericCoefficient(int order, const Partition& partition, int sampleSize)
{
    const QMap<int, int> counts = partition.multiplicities();
    const int blocks = partition.blockCount();
    const long double factor = static_cast<long double>(multinomialBlockFactor(order, counts));

    long double falling = 1.0L;
    for (int i = 1; i <= blocks - 1; ++i) {
        falling *= static_cast<long double>(sampleSize - i);
    }

    long double bracket = 0.0L;
    for (auto it = counts.cbegin(); it != counts.cend(); ++it) {
        const int q = it.key();
        const int count = it.value();
        const long double sign = ((order - q) % 2 == 0) ? 1.0L : -1.0L;
        bracket += sign * static_cast<long double>(count) * std::pow(sampleSize - 1.0L, q);
    }
    bracket += (order % 2 == 0 ? 1.0L : -1.0L) * static_cast<long double>(sampleSize - blocks);

    return factor * falling * bracket / std::pow(static_cast<long double>(sampleSize), order);
}

QString MomentEngine::formatNumber(long double value)
{
    if (std::fabs(value) < 1e-14L) {
        value = 0.0L;
    }
    return QString::number(static_cast<double>(value), 'g', 10);
}
