#include "gdarap.h"

#include <Mesh.hpp>
#include <ArrayMesh.hpp>

#include <igl/remove_unreferenced.h>
#include <igl/point_mesh_squared_distance.h>

#include <map>

using namespace godot;

void GDArap::_register_methods()
{
    register_method("_process", &GDArap::_process);
    register_method("_ready", &GDArap::_ready);
    
    register_property("pinned_indices", &GDArap::pinned_indices, PoolIntArray());
    register_property("max_iterations", &GDArap::max_iter, 10);
    register_property("energy", &GDArap::energy, 3);
    register_property("youngs_modulus", &GDArap::youngs_modulus, 1.0);
    register_property("init_randomness", &GDArap::randomness, 0.0);
    register_property("with_dynamics", &GDArap::with_dynamics, true);
    register_property("gravity", &GDArap::gravity, 9.8);
    register_property("coeff_of_restitution", &GDArap::CoR, 1.1);
    register_property("ground_height", &GDArap::ground_height, -0.8);

    register_method("_on_MI_value_changed", &GDArap::set_max_iter);
    register_method("_on_E_value_changed", &GDArap::set_energy);
    register_method("_on_YM_value_changed", &GDArap::set_youngs_modulus);
    register_method("_on_G_value_changed", &GDArap::set_gravity);
    register_method("_on_CoR_value_changed", &GDArap::set_CoR);
    register_method("_on_GH_value_changed", &GDArap::set_ground_height);
}

GDArap::GDArap(){}

GDArap::~GDArap(){}

void GDArap::_init()
{
    // initialize values
    is_initialized = false;
    max_iter = 10;
    energy = 3;
    youngs_modulus = 1;
    randomness = 0.0;
    with_dynamics = true;
    gravity = 9.8;
    CoR = 1.1;
    ground_height = -0.8;
}

void GDArap::_ready()
{
    is_initialized = false;

    Ref<Mesh> mesh = get_mesh();
    if (mesh.is_null())
    {
        return;
    }

    // obtain V & F from mesh, only handle surface 0 here

    // obtain verts
    mesh_arrays = mesh->surface_get_arrays(0);
    PoolVector3Array verts = mesh_arrays[Mesh::ARRAY_VERTEX];
    Eigen::MatrixXd V0(verts.size(), 3);
    for (unsigned int vi = 0; vi < V0.rows(); ++vi)
    {
        V0.row(vi) << verts[vi][0], verts[vi][1], verts[vi][2];
    }

    // verts appear to be a non-indexed array i.e. shared vertices are duplicated in the array, see
    // https://docs.godotengine.org/en/stable/tutorials/content/procedural_geometry/index.html#surface-array
    // mapping indices I for unique vertices based on distances.
    {
        Eigen::VectorXi E = Eigen::VectorXi::LinSpaced(V0.rows(), 0, V0.rows() - 1);
        Eigen::VectorXd sqrD;
        Eigen::MatrixXd C;
        igl::point_mesh_squared_distance(V0, V0, E, sqrD, I, C);
    }

    // obtain faces for F
    PoolIntArray faces = mesh_arrays[Mesh::ARRAY_INDEX];
    Eigen::MatrixXi F(faces.size() / 3, 3);
    for (unsigned int fi = 0; fi < F.rows(); ++fi)
    {
        unsigned int ff = fi * 3;
        // pointing vertices to V0 instead of verts
        F.row(fi) << I(faces[ff]), I(faces[ff+1]), I(faces[ff+2]);
    }

    { // remove unreferenced vertices in V0 & re-map I to the cleaned-up V
        Eigen::VectorXi I_, J_;
        igl::remove_unreferenced(V0, Eigen::MatrixXi(F), V, F, I_, J_);
        std::map<int, int> idx_map;
        for (int j = 0; j < J_.size(); ++j)
        {
            idx_map[ J_(j) ] = j;
        }
        for (int i = 0; i < I.size(); ++i)
        {
            I(i) = idx_map[ I(i) ];
        }
    }

    // set pinned indices as ARAP boundary conditions
    const size_t n = V.rows();
    const unsigned int n_pins = pinned_indices.size();
    Eigen::VectorXi b(n_pins);
    for (unsigned int i = 0; i < n_pins; ++i)
    {
        b(i) = std::max<int>(0, std::min<int>(n - 1, pinned_indices[i]));
    }

    // precompute ARAP
    arap_data.with_dynamics = with_dynamics;
    arap_data.max_iter = max_iter;
    arap_data.energy = static_cast<igl::ARAPEnergyType>(energy);
    arap_data.ym = 0.01 * youngs_modulus;
    if(!igl::arap_precomputation(V, F, V.cols(), b, arap_data))
    {
        return;
    }

    // init ARAP physics conditions
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);
    arap_data.f_ext =  M * Eigen::RowVector3d(0, -gravity, 0).replicate(n, 1);
    // adding random initial velocities
    arap_data.vel = randomness * Eigen::MatrixXd::Random(n, 3);

    is_initialized = true;
}

void GDArap::_process(float delta)
{
    if (!is_initialized)
    {
        return;
    }

    // setting boundary condition from pinned indices
    const size_t n = V.rows();
    const unsigned int n_pins = pinned_indices.size();
    Eigen::MatrixXd bc(n_pins, 3);
    for (unsigned int i = 0; i < n_pins; ++i)
    {
        const int bi = std::max<int>(0, std::min<int>(n - 1, pinned_indices[i]));
        bc.row(i) << V.row(bi);
    }

    // ARAP solve step
    arap_data.h = delta;
    arap_data.max_iter = max_iter;
    arap_data.energy = static_cast<igl::ARAPEnergyType>(energy);
    arap_data.ym = 0.01 * youngs_modulus;
    arap_data.f_ext =  M * Eigen::RowVector3d(0, -gravity, 0).replicate(n, 1);
    igl::arap_solve(bc, arap_data, V);

    // collision with ground
    for (unsigned int vi = 0; vi < V.rows(); ++vi)
    {
        const int y = 1;
        if(V(vi, y) < ground_height)
        {
            V(vi, y) = ground_height - V(vi, y) + ground_height;
            arap_data.vel(vi, y) = - arap_data.vel(vi, y) / CoR;
        }
    }

    // update mesh with V
    PoolVector3Array verts = mesh_arrays[Mesh::ARRAY_VERTEX];
    for (unsigned int ii = 0; ii < verts.size(); ++ii)
    {
        const auto& v = V.row(I(ii));
        verts.set(ii, Vector3(v(0), v(1), v(2)));
    }
    mesh_arrays[Mesh::ARRAY_VERTEX] = verts;

    ArrayMesh* array_mesh = ArrayMesh::_new();
    array_mesh->add_surface_from_arrays(Mesh::PRIMITIVE_TRIANGLES, mesh_arrays);
    set_mesh(array_mesh);
}