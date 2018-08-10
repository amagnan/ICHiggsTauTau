#include "TFile.h"
#include <iostream>
#include <fstream>
using namespace std;
int againstZombie_data(){

  fstream fList("list_data");

  std::ofstream outfile;
  outfile.open("failed_files_DATA_reminiaod_MuEGClean.dat");
  int i=0;
  char buf[10000];

  string prefix("root://gfe02.grid.hep.ph.ic.ac.uk:1095//store/user/rdimaria/170201_DATA");
  while( !fList.eof() ){
    fList.getline(buf, 10000);
    cout << i++ << " " << buf << endl;
    fstream fTreeList(buf);
    char bufTree[10000];
    while( !fTreeList.eof() ){
      fTreeList.getline(bufTree, 10000);
      string sBuf(bufTree);
      string name = prefix+bufTree; 
      cout << name << endl;

      TFile *f = TFile::Open(name.c_str());
      if (!f) outfile << bufTree << std::endl;
    }
  }
  outfile.close();
  return 0;
}
