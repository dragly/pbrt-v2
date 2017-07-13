#include "Sprite.h"
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <fstream>
#include <math.h>

Sprite::Sprite(const std::string& dataDirectory, const std::string& header,
               const uint &dim) {
    spriteDirectory = dataDirectory;
    headerFile = header;
    gridDim = dim;

    xmin = std::numeric_limits<float>::max();
    ymin = std::numeric_limits<float>::max();
    zmin = std::numeric_limits<float>::max();
    xmax = std::numeric_limits<float>::min();
    ymax = std::numeric_limits<float>::min();
    zmax = std::numeric_limits<float>::min();

    // Build the event sprite
    BuildFluorescenceMap(spriteDirectory, headerFile);
    ComputeBoundingBox();
}

void Sprite::ReadHeaderData(const std::string &header) {
    std::ifstream file(header.c_str());
    std::string line;
    if (file.is_open()) {
        while (std::getline(file, line))
        {
            ParseIntParameter(line, "EventsCount", eventsCount);
            ParseFloatParameter(line, "XCenter", xCenter);
            ParseFloatParameter(line, "YCenter", yCenter);
            ParseFloatParameter(line, "ZCenter", zCenter);
            ParseFloatParameter(line, "XCOI", xCOI);
            ParseFloatParameter(line, "YCOI", yCOI);
            ParseFloatParameter(line, "ZCOI", zCOI);
            ParseFloatParameter(line, "AABBWidth", width);
            ParseFloatParameter(line, "AABBHeight", height);
            ParseFloatParameter(line, "AABBDepth", depth);
            ParseStringParameter(line, "VSDPositionFile", positionFile);
            ParseStringParameter(line, "VSDIntensityFile", intensityFile);
            ParseFloatParameter(line, "TimeStep", timeStep);
        }
        file.close();
    } else {
        printf("The specified header [%s] can't be opened \n", header.c_str());
    }

    printf("Sprite: [%s] \n", header.c_str());
    printf("\t EventsCount %u \n", eventsCount);
    printf("\t Center [%f %f %f] \n", xCenter, yCenter, zCenter);
    printf("\t Center of Interest [%f %f %f] \n", xCOI, yCOI, zCOI);
    printf("\t Bounding Box [%.5f %.5f %.5f] \n", width, height, depth);
    printf("\t VSDPositionFile [%s] \n", positionFile.c_str());
    printf("\t VSDIntensityFile [%s] \n", intensityFile.c_str());
}

void Sprite::BuildFluorescenceMap(const std::string &dataDirectory,
                                  const std::string &pshFile) {

    // Read the header file to extract the data
    const std::string pshFilePath = dataDirectory + "/" + pshFile;
    printf("Path %s \n", pshFilePath.c_str());
    ReadHeaderData(pshFilePath);

    // Open the data file and read the fluorescence events
    std::ifstream positionFileStream, intensityFileStream;
    const std::string pspFilePath = dataDirectory + "/" + positionFile;
    const std::string psiFilePath = dataDirectory + "/" + intensityFile;
    printf("Sources [%s][%s] \n", pspFilePath.c_str(), psiFilePath.c_str());
    positionFileStream.open(pspFilePath.c_str(), std::ios::binary);
    intensityFileStream.open(psiFilePath.c_str(), std::ios::binary);
    if(positionFileStream.is_open() && intensityFileStream.is_open())
    {
        for(uint i = 0; i < eventsCount; i++)
        {
            float x, y, z, intensity;
            positionFileStream.read(reinterpret_cast<char*>(&x), sizeof(x));
            positionFileStream.read(reinterpret_cast<char*>(&y), sizeof(y));
            positionFileStream.read(reinterpret_cast<char*>(&z), sizeof(z));
            intensityFileStream.read(reinterpret_cast<char*>(&intensity),
                                     sizeof(intensity));
            fluorescenceSources.push_back(FluorescenceEvent(
                                              Point(x - xCOI,
                                                    y - yCOI,
                                                    z - zCOI),
                                              intensity));
        }

        positionFileStream.close();
        intensityFileStream.close();
    }
    printf("The source data has %d fluorescent events \n",
           uint(fluorescenceSources.size()));
}

void Sprite::ComputeBoundingBox() {
    for(uint i = 0; i < eventsCount; i++) {
        FluorescenceEvent& event = fluorescenceSources[i];
        if(event.point.x > xmax) xmax = event.point.x;
        if(event.point.y > ymax) ymax = event.point.y;
        if(event.point.z > zmax) zmax = event.point.z;

        if(event.point.x < xmin) xmin = event.point.x;
        if(event.point.y < ymin) ymin = event.point.y;
        if(event.point.z < zmin) zmin = event.point.z;
    }
}

void Sprite::WriteSprite(std::string outputDir) {
    size_t lastindex = headerFile.find_last_of(".");
    std::string fileName = headerFile.substr(0, lastindex);
    std::string psFile = outputDir + "/" + fileName + ".ps";
    std::ofstream psStream(psFile.c_str(), std::ios::out);
    for(uint i = 0; i < eventsCount; i++)
    {
        FluorescenceEvent& event = fluorescenceSources[i];
        psStream.write((char *)&event.point.x, 4);
        psStream.write((char *)&event.point.y, 4);
        psStream.write((char *)&event.point.z, 4);
        psStream.write((char *)&event.intensity, 4);

    }
    psStream.close();
}

void Sprite::SampleSpriteToVolume(std::string outputDir) {
    printf("Volume: \n");
    printf("\t Extent [%f %f %f] \n", width, height, depth);
    printf("\t pMin [%f %f %f] \n", xmin, ymin, zmin);
    printf("\t pMax [%f %f %f] \n", xmax, ymax, zmax);

    float deltaX = xmax - xmin;
    float deltaY = ymax - ymin;
    float deltaZ = zmax - zmin;
    float nwidth, nheight, ndepth;
    float largestDim = width;
    if(height > largestDim) largestDim = height;
    if(depth > largestDim) largestDim = depth;
    nwidth = width / largestDim;
    nheight = height / largestDim;
    ndepth = depth / largestDim;
    printf("\t Normalized Extent [%f %f %f] \n", nwidth, nheight, ndepth);

    uint nx, ny, nz;
    nx = (uint) (nwidth * (gridDim));
    ny = (uint) (nheight * (gridDim));
    nz = (uint) (ndepth * (gridDim));
    uint nxyz = nx * ny * nz;
    printf("\t Grid Size [%u %u %u] \n", nx, ny, nz);

    // Create the volume grid, initialize it and sample the events to it
    volume = new float[nxyz];
    for(int i = 0; i < nxyz; i++) volume[i] = 0;
    for(uint i = 0; i < eventsCount; i++) {
        FluorescenceEvent& event = fluorescenceSources[i];
        uint vx = (uint) (nx * (floor(event.point.x - xmin) / deltaX));
        uint vy = (uint) (ny * (floor(event.point.y - ymin) / deltaY));
        uint vz = (uint) (nz * (floor(event.point.z - zmin) / deltaZ));
        uint idx = vx + nx * vy + nx * ny * vz;

        // TODO: fix with the actual values
        float value = 0;
        if (event.intensity < 0)
            value = -1.f * event.intensity;
        else
            value = event.intensity;
        volume[idx] += value ;
    }

    // Normalize the grid
    float maxValue = 0.f;
    for(uint i = 0; i < nxyz; i++)
        if(volume[i] > maxValue) maxValue = volume[i];
    // Set the maximum value to 255
    for(uint i = 0; i < nxyz; i++) {
        // TODO: Fix with the actual values
        if(volume[i] > 0) {
            volume[i] =  255.0 * ((volume[i] / maxValue));
        }
    }

    // Write the volume (.hdr/.img or .raw)
    size_t lastindex = headerFile.find_last_of(".");
    std::string fileName = headerFile.substr(0, lastindex);
    std::string hdr = outputDir + "/" + fileName + ".hdr";
    std::string img = outputDir + "/" + fileName + ".img";

    std::ofstream hdrStream(hdr.c_str(), std::ios::out);
    hdrStream << nx << " " << ny << " " << nz << "\n";
    hdrStream.close();

    std::ofstream imgStream(img.c_str(), std::ios::out | std::ios::binary);
    for(int i = 0; i < nxyz; i++) {
        int value = int(volume[i]);
        imgStream.write((char *)&value, 1);
    }
    imgStream.close();

    delete volume;
}


void Sprite::ParseStringParameter(const std::string& line,
                                  const std::string& parameter,
                                  std::string& value) {
    if (line.find(std::string(parameter).append("=")) != std::string::npos) {
        std::istringstream iss(line);
        std::string token;
        int counter = 0;
        while(std::getline(iss, token, '=')) {
            if (counter++ == 1) {
                value = token;
            }
        }
    }
}

void Sprite::ParseIntParameter(const std::string& line,
                               const std::string& parameter,
                               uint& value) {
    if (line.find(std::string(parameter).append("=")) != std::string::npos) {
        std::istringstream iss(line);
        std::string token;
        int counter = 0;
        while(std::getline(iss, token, '=')) {
            if (counter++ == 1) {
                value = ::atoi(token.c_str());
            }
        }
    }
}

void Sprite::ParseFloatParameter(const std::string& line,
                                 const std::string& parameter,
                                 float& value) {
    if (line.find(std::string(parameter).append("=")) != std::string::npos) {
        std::istringstream iss(line);
        std::string token;
        int counter = 0;
        while(std::getline(iss, token, '=')) {
            if (counter++ == 1) {
                value = ::atof(token.c_str());
            }
        }
    }
}
