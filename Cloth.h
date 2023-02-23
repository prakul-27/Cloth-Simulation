#include "Particle.h" 
#include "Spring.h"
#include <bits/stdc++.h>

#define ACCELERATION_OF_GRAVITY 9.80f
#define MESH_SIZE 15 

class Cloth{
public:
    // info about mesh particle (i, j)
    std::vector<std::vector<Particle>>particles;
    // info about all springs connected from/to (i, j) [vertex set]
    std::vector<std::vector<std::vector<Spring>>> springs;         
    // info about springs [edge set]
    std::vector<Spring>structural_springs;
    std::vector<Spring>shear_springs;
    std::vector<Spring>flexion_springs;
    

    // set initial state of the cloth 
    Cloth();



    // simulate the Cloth movement
    void Simulate();
    // compute net gravitational force on mesh particle (i, j)
    Eigen::Vector3d computeNetGravitationalForce(int i, int j);
    // compute net internal force due to springs on mesh particle (i, j)
    Eigen::Vector3d computeNetInternalForce(int i, int j);
    
    // accumulate net internal forces
    void accumulateNetInternalForces();
    // compute net damping force on mesh particle (i, j)
    Eigen::Vector3d computeNetDampingForce(int i, int j);
    // compute net viscous force on mesh particle (i, j)
    Eigen::Vector3d computeNetViscousForce(int i, int j);
    
    // compute unit normal on the surface at point (i,j)
    Eigen::Vector3d computeNormal(int i, int j);
    void computeNormal();    



    // render the mesh onto the screen
    void Draw();

    //-----new added--------
   // void ApplyEulerMethod(void);
	//void ApplyDynamicInverseProcedure(void);

    // debugging functions
    void printState(int time);
};