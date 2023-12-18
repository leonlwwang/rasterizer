#ifndef RASTERIZER_H
#define RASTERIZER_H

#include <Eigen/Dense>
#include <optional>
#include <fstream>

using std::string, Eigen::VectorXd;

/**
 * Draws a triangle.
 * @param first the initial index position
 * @param count a multiple of 3
*/
void drawArraysTriangles(unsigned first, unsigned count);

/**
 * Divides RGB values by w if hyperbolic interpolation is on
 * @param vec the color vector
 * @param w the w value from position vector
*/
VectorXd divideRGB(VectorXd vec, double w);

/**
 * Transforms the given vector by division of w and viewport
 * @param vec the position vector
 * @returns the transformed vector
*/
VectorXd viewportAndDivision(VectorXd vec);

/**
 * Draws the given vector onto the image representing a pixel
 * @param v the vector representing a pixel
*/
void draw(VectorXd &v);

/**
 * DDA algorithm to find all whole points on a line segment.
 * @param a vector
 * @param b vector
 * @param axis which axis to perform algorithm on (0 -> x, 1 -> y)
*/
std::optional<std::pair<VectorXd, VectorXd>> dda(VectorXd a, VectorXd b, bool axis);

/**
 * Performs the scanline algorithm to find all points inside a triangle.
 * @param P triangle vertex (8D vector)
 * @param Q triangle vertex (8D vector)
 * @param R triangle vertex (8D vector)
 * @param image the image to draw on
*/
void scanLine(const VectorXd &P, const VectorXd &Q, const VectorXd &R);

/**
 * Goes through the image data, processes any mode/state changes, loads
 * any objects found in the data and saves each of them in a vector of objects
 * @param file .txt file with image data
*/
void processCommands(std::ifstream& file, bool debug = false);

/**
 * Initializes global variables for the image data like width, height, and file name
 * @param file .txt file with image data
*/
void setImageGlobals(std::ifstream& file, bool debug = false);

/**
 * Main function that generates the image via raytracing, used in main.cpp
 * @param file .txt file with image data
 * @returns a pair of pairs, first pair is (scene data, scene name), second pair is (x, y) of scene dimensions
*/
std::pair<std::pair<std::vector<unsigned char>, string>, std::pair<int, int>> generateImage(std::ifstream& file);

/**
 * Utility function that opens the image data .txt file
 * @param argc command line argument
 * @param argv command line argument
 * @param debug prints the file details
*/
std::ifstream openFile(int argc, char** argv, bool debug = false);

#endif