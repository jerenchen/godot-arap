#ifndef GDARAP_H
#define GDARAP_H

#include <Godot.hpp>
#include <MeshInstance.hpp>

#include <Eigen/Dense>
#include <igl/arap.h>


namespace godot {

class GDArap : public MeshInstance {
    GODOT_CLASS(GDArap, MeshInstance)

private:
    bool is_initialized;
    Array mesh_arrays;

    // ARAP data/parameters
    igl::ARAPData arap_data;
    int max_iter;
    int energy;
    double youngs_modulus;
    double randomness;
    bool with_dynamics;
    double gravity;
    double CoR; // Coefficient of Restitution
    double ground_height;
    PoolIntArray pinned_indices;

    Eigen::SparseMatrix<double> M; // mass matrix
    Eigen::VectorXi I; // indices mapping verts (non-indexed) to V
    Eigen::MatrixXd V; // unique set of vertices

public:
    // Godot node funcs
    static void _register_methods();

    GDArap();
    ~GDArap();

    void _init();
    void _ready();
    void _process(float delta);

    // setter funcs for receiving signals
    void set_max_iter(int value)
    {
	   max_iter = value;
    }
    void set_energy(int value)
    {
	   energy = std::max<int>(0, std::min<int>(value, 3));
    }
    void set_youngs_modulus(float value)
    {
	   youngs_modulus = (double)value;
    }
    void set_gravity(float value)
    {
	   gravity = (double)value;
    }
    void set_CoR(float value)
    {
	   CoR = (double)value;
    }
    void set_ground_height(float value)
    {
	   ground_height = (double)value;
    }
};

}

#endif