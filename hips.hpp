#pragma once

#include <QString>
#include <QImage>

constexpr int HIPS_TILE_SIZE = 512;
inline const QString hipsInitialExt = "tiff"; // Should be uncompressed for saving and reloading speed
QString formatDate();
QString md5sum(QString const& filePath);
void createLowerOrderTile(const int order, const int pix, const QString& outDir);
void generateLowerOrderTiles(const int orderMax, const QString& outDir, bool savePartialTiles = false); // Merge the deepest-level tiles to create the ones with lower detail level
void convertTiles(const QString& finalExt, const QString& formatName, const int orderMax, const QString& outDir);
void hipsSaveProperties(QString const& outDir, const int orderMax, QString const& imgFormat, QString const& surveyTitle,
                        QString const& surveyType, QString const& description, QString const& frame, QString const& obs_copyright,
                        QString const& hips_copyright, QString const& creator, QString const& hipsStatus);
