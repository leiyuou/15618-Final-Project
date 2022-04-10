/**
 * 15618 Project shortest.h
 * 
 * Xinqi Wang, Yuou Lei
 */

#include <omp.h>
#include <bits/stdc++.h>

typedef std::pair<int,int> adjPair;
typedef std::pair<int,int> distPair;

 struct compare  
 {  
   bool operator()(const distPair& d1, const distPair& d2)  
   {  
       return d1.second > d2.second;  
   }  
 }; 