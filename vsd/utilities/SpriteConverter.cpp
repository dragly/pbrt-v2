#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <unistd.h>
#include <ctype.h>
#include <getopt.h>
#include "Sprite.h"

int
main(int argc, char** argv) {

    std::string inputDir;
    std::string outputDir;
    std::string pshFile;
    int volDim = 0;
    int c;
    while ((c = getopt (argc, argv, "hi:o:p:d:")) != -1)
        switch (c)
        {
        case 'h':
            printf("arguments: -i <input directory> -p <psh file> -o "
                   "<output directory> -d <grid dimensions> \n");
            return EXIT_SUCCESS;
        case 'i':
            inputDir = std::string(optarg);
            break;
        case 'o':
            outputDir = std::string(optarg);
            break;
        case 'p':
           pshFile = std::string(optarg);
            break;
        case 'd':
            volDim = atoi(optarg);
            break;
        default:
            abort ();
        }

    printf("Input: %s \nOutput: %s \npsh file: %s \nGrid dim: %d \n",
           inputDir.c_str(), outputDir.c_str(), pshFile.c_str(), volDim);

    // Read the sprite file
    Sprite* sprite = new Sprite(inputDir, pshFile, volDim);
    sprite->SampleSpriteToVolume(outputDir);

    // Double check the results
    sprite->WriteSprite(outputDir);

    // Cleaning
    delete sprite;

    return EXIT_SUCCESS;
}
