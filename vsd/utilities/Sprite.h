#ifndef SPRITE_H
#define SPRITE_H

#include <string>
#include <vector>
#include "Common.h"

struct Point {
    float x;
    float y;
    float z;
    Point(const float& xx, const float& yy, const float& zz)
        : x(xx), y(yy), z(zz) { }
};

class FluorescenceEvent {
public:
    FluorescenceEvent (const Point p, const float i):
        point(p), intensity(i) { }
    Point point;
    float intensity;
};


class Sprite
{
public:
    Sprite(const std::string& dataDirectory, const std::string& header,
           const uint& dim);
    void ReadHeaderData(const std::string& header);
    void ComputeBoundingBox();
    void SampleSpriteToVolume(std::string outputDir);
    void WriteSprite(std::string outputDir);
    void BuildFluorescenceMap(const std::string& dataDirectory,
                              const std::string& pshFile);
    void ParseIntParameter(const std::string& line,
                           const std::string& parameter,
                           uint& value);
    void ParseFloatParameter(const std::string& line,
                             const std::string& parameter,
                             float& value);
    void ParseStringParameter(const std::string& line,
                              const std::string& parameter,
                              std::string& value);
private:
    uint eventsCount;
    float xCenter, yCenter, zCenter;
    float xCOI, yCOI, zCOI;
    float width, height, depth;
    float xmin, ymin, zmin;
    float xmax, ymax, zmax;
    uint gridDim;
    std::string spriteDirectory;
    std::string headerFile;
    std::string positionFile;
    std::string intensityFile;
    float timeStep;
    std::vector <FluorescenceEvent> fluorescenceSources;
    float* volume;
};

#endif // SPRITE_H
