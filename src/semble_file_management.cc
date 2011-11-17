// semble_file_management.cc -
//
// Saturday, October 22 2011
//

#include"semble_file_management.h"


namespace SEMBLE
{

  namespace SEMBLEIO
  {
    std::string getPath(void)
    {
      char cCurrentPath[FILENAME_MAX];

      if(!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
        {
          std::cout << "Couldnt get path in " << __PRETTY_FUNCTION__ << __FILE__ << __LINE__ << std::endl;
          exit(1);
        }

      cCurrentPath[sizeof(cCurrentPath) - 1] = '\0'; //stick in a null character
      std::string dum(cCurrentPath);
      return  dum += std::string("/");
    }


    void makeDirectoryPath(const std::string &s)
    {
      if(!!!(access(s.c_str(), 0) == 0))  //check if it already exists, theres no reason to use extra system cmds if we dont have to since its bad
        {                                 //NB, access returns true if there was a file OR directory with that name, would be nice to check if its a directory
          std::cout << "Making path:" << s << std::endl;
          std::string cmd = "mkdir -p ";
          cmd += s;
          system(cmd.c_str());
        }
    }

  }
}
