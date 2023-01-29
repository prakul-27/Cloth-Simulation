#include "Cloth.h"
#include <unistd.h> 

int main() {    
    Cloth cloth;
    int t = 0;
    while(1) {
        cloth.Simulate();            
        cloth.printState(t);
        // return 0;
        t++;
        // if(t == 2) break;
    }       
}