#pragma once
#include "global.h"
#include <fstream>
#include <cstdio>
#include <iostream>
#include <highfive/highfive.hpp>

bool saveVTK_Scalar(matrix3d<double> data,
    int t,
	long int xS, long int yS, long int zS,
	long int xE, long int yE, long int zE,
	double dx, double dy, double dz,
    const std::string& fileRoot, const std::string& varName,
    int appendFlag);

bool saveVTK_Vector(matrix3d<double> data1, matrix3d<double> data2, matrix3d<double> data3,
    int t,
	long int xS, long int yS, long int zS,
	long int xE, long int yE, long int zE,
	double dx, double dy, double dz,
    const std::string& fileRoot, const std::string& varName,
    int appendFlag);