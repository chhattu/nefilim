// Copyright (C)  2012 Zhang Initative Research Unit
// SegmentProfile.cc -description
// written by Rojan Shrestha 

#include<QueryPFImpl.hh>

#include<stdlib.h>
#include<iostream>
#include<string>
#include<map>

using namespace std;

// what could be error arisen situations?
void option_parser(int argc, char* argv[], map<string, string> &args) {
  int i = 1;
  while(i < argc) { 
    args[argv[i]] = argv[i+1]; 
    i+=2;
  }
}

void usage() {
  cout<<"Usage:"<<endl;
  cout<<"  QueryProfileVector -qv path2queryvector -qc  \
            path2querycoord -dv path2dbvector -dc path2dbcoord"<<endl;
  cout<<"   -qv path to query vector file."<<endl;
  cout<<"   -qc path to query coordinates file."<<endl;
  cout<<"   -dc path to database vector file."<<endl;
  cout<<"   -dv path to database coordinate file."<<endl;
  exit(1); 
}

int main(int argc, char* argv[]) {
  map<string, string> args;  
  option_parser(argc, argv, args);
  if(args.empty()) usage();

  if(args["-qc"].empty() || args["-qv"].empty() || 
     args["-dc"].empty() || args["-dv"].empty())  {
    usage();
  }

  QueryPFImpl query_pf_impl; 
  if(query_pf_impl.set_query_vector(args["-qv"], args["-qc"]) && 
     query_pf_impl.set_database_vector(args["-dv"], args["-dc"])) {
    query_pf_impl.query();
  }
}


