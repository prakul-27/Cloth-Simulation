#include "Particle.h"

// specify position as (x, y, z)
Particle::Particle(double x, double y, double z) {
    position = Eigen::Vector3d(x, y, z); 
    velocity.setZero();
    normal.setZero(); 
    netForce.setZero();
}

// specify position as vector
Particle::Particle(Eigen::Vector3d position) {
    this->position = position;
    velocity.setZero();
    normal.setZero();
    netForce.setZero();
}

// set everything to origin
Particle::Particle() {
    position.setZero();
    velocity.setZero();
    normal.setZero();
    netForce.setZero();
}

// render on screen
void Particle::Draw() {

}

// zero the forces acting on it
void Particle::clearForces() {
    netForce.setZero();
}