#pragma once

#include <QMap>
#include <QString>
#include <QVector>

#include <optional>

struct FormulaSection
{
    QString title;
    QString latex;
};

struct Partition
{
    QVector<int> parts;

    QMap<int, int> multiplicities() const;
    int blockCount() const;
    QString labelLatex() const;
    QString productLatex() const;
};

struct MomentTerm
{
    Partition partition;
    QString coefficientSymbolLatex;
    QString coefficientLatex;
    QString productLatex;
    QString symbolicTermLatex;
    QString numericTermLatex;
    std::optional<long double> numericCoefficient;
};

class MomentEngine
{
public:
    static QVector<Partition> partitionsWithoutOnes(int order);
    static QVector<MomentTerm> terms(int order, std::optional<int> sampleSize);
    static QVector<FormulaSection> formulaSections(int order, const QVector<MomentTerm>& terms, std::optional<int> sampleSize);
    static QString fullLatex(int order, const QVector<MomentTerm>& terms, std::optional<int> sampleSize);

private:
    static void generatePartitions(int remaining, int maxPart, QVector<int>& current, QVector<Partition>& out);
    static qint64 factorial(int value);
    static qint64 multinomialBlockFactor(int order, const QMap<int, int>& multiplicities);
    static QString coefficientLatex(int order, const Partition& partition);
    static QString bracketLatex(int order, const QMap<int, int>& multiplicities, int blockCount);
    static QString fallingFactorLatex(int start, int count);
    static QString definitionsLatex(int order);
    static QString forwardLatex(int order, const QVector<MomentTerm>& terms, std::optional<int> sampleSize);
    static QString inverseLatex(int order, const QVector<MomentTerm>& terms);
    static QString applicationsLatex(int order);
    static QString estimatorProductLatex(const Partition& partition);
    static QString joinTerms(const QStringList& terms);
    static const MomentTerm* leadingTerm(const QVector<MomentTerm>& terms, int order);
    static long double numericCoefficient(int order, const Partition& partition, int sampleSize);
    static QString formatNumber(long double value);
};
