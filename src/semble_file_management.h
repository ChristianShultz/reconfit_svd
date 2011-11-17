#ifndef SEMBLE_FILE_MANAGEMENT_H_H_GUARD
#define SEMBLE_FILE_MANAGEMENT_H_H_GUARD

#include<iostream>
#include <stdio.h>  /* defines FILENAME_MAX -- 1024 usually I think*/
#include<stdlib.h>
#include<string>

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

/*
  Two file management functions to find the directory where a program is being run and
  another to create a path to stash files into
*/

namespace SEMBLE
{
  namespace SEMBLEIO
  {
    //return a path to where I am (pwd)
    std::string getPath(void);

    //try to make the path s, it isn't recursive so you should make each folder from top to bottom if you want multiple levels
    void makeDirectoryPath(const std::string &s);

  }
}

#endif
