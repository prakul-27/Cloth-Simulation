#include "Spring.h"

Spring::Spring(int type, Particle *p1, Particle *p2) {
    this->type = type;    
    this->p1 = p1; 
    this->p2 = p2;
}

Spring::Spring(int type) {
    this->type = type;
}

Eigen::Vector3d Spring::computeSpringForce() {    
    float spring_natural_length = type == TYPE_STRUCTURAL ? STRUCTURAL_NATURAL_LENGTH:
                        type == TYPE_SHEAR ? SHEAR_NARTURAL_LENGTH:FLEXION_NATURAL_LENGTH;
    float spring_consant = type == TYPE_STRUCTURAL ? STRUCTURAL_STIFFNESS:
                        type == TYPE_SHEAR ? SHEAR_STIFFNESS:FLEXION_STIFFNESS;        
    Eigen::Vector3d particle_displacement = p2->position - p1->position;
    Eigen::Vector3d normal = particle_displacement/particle_displacement.norm();
    Eigen::Vector3d spring_displacement = spring_natural_length*normal;
    Eigen::Vector3d diff = particle_displacement - spring_displacement;    
    return diff*spring_consant;
}