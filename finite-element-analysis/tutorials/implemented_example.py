from finiteelementanalysis import pre_process as pre
from finiteelementanalysis import pre_process_demo_helper_fcns as pre_demo
from finiteelementanalysis.solver import hyperelastic_solver
from finiteelementanalysis import visualize as viz
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


# for saving files later
tutorials_dir = Path(__file__).parent

# --- Beam geometry ---
L = 60.0   # length in x
H = 1  # height in y
nx = 60 #40    # number of elements along length
ny = 2     # number of elements along height

ele_type = "D2_nn8_quad"  # 2D, 4-node quadrilateral (linear)
ndof = 2                  # 2 DOFs per node (x, y)

# Generate a rectangular mesh
coords, connect = pre.generate_rect_mesh_2d(ele_type, 0.0, 0.0, L, H, nx, ny)
# coords: shape (n_nodes, 2)
# connect: shape (n_nodes_per_elem, n_elems)

# --- Identify boundaries ---
boundary_nodes, boundary_edges = pre.identify_rect_boundaries(
    coords, connect, ele_type, x_lower=0.0, x_upper=L, y_lower=0.0, y_upper=H
)

# 1) Clamp the left edge: fix x- and y-displacements = 0
fixed_left = pre.assign_fixed_nodes_rect(boundary_nodes, "left", 0.0, 0.0)
fixed_right = pre.assign_fixed_nodes_rect(boundary_nodes, "right", L, 0.0)


# 2) Uniform downward traction on the top edge (y=H)
# Let q be negative in the y-direction
q = -0.8  # load per unit length in x
# For a 2D plane strain problem, this is a traction (tx, ty) = (0, q)
dload_info = pre.assign_uniform_load_rect(boundary_edges, "top", 0.0, q)

# Combine boundary conditions
fixed_nodes = np.hstack((fixed_left, fixed_right))


# --- Material properties ---
E = 100000.0
nu = 0.3
mu = E / (2.0 * (1.0 + nu))
kappa = E / (3.0 * (1.0 - 2.0 * nu))


material_props = np.array([mu, kappa])
print(f"Material properties: mu={mu:.3f}, kappa={kappa:.3f}")

# Number of incremental load steps
nr_num_steps = 1

# --- Solve with your hyperelastic solver ---
displacements_all, nr_info_all = hyperelastic_solver(
    material_props,
    ele_type,
    coords.T,      # shape (2, n_nodes)
    connect.T,     # shape (n_nodes_per_elem, n_elems)
    fixed_nodes,
    dload_info,
    nr_print=True,
    nr_num_steps=nr_num_steps,
    nr_tol=1e-10,
    nr_maxit=300,
)

final_disp = displacements_all[-1]  # shape: (n_nodes*ndof,)

# --- Compute the tip displacement from the FEA result ---
# We'll pick a node near x=L, y=H/2
tip_node = None
tol = 1e-2
for i, (x, y) in enumerate(coords):
    if abs(x - 0.5*L) < tol and abs(y - H/2) < H/(2*ny):
        tip_node = i
        break
if tip_node is None:
    raise ValueError("Could not find tip node near x=L, y=H/2.")

mid_disp_y = final_disp[ndof*tip_node + 1]  # the y-component of displacement

# --- Compare with Euler–Bernoulli formula for small deflection ---
mu = material_props[0]
# We'll guess E = 3 mu (plane stress, near incompressible) or 2 mu(1 + nu) for plane strain, etc.
# --- Analytical Solution ---
# For a cantilever beam under uniformly distributed load q:
# Euler–Bernoulli tip deflection: w(L) = q * L^4 / (8 * E * I)
# E_eff = E * (1 - nu) / ((1 + nu) * (1 - 2*nu))
E_eff = E / (1 - nu ** 2.0)
I = H ** 3 / 12.0
w_analytical = q * L ** 4 / (8.0 * E_eff * I)

print(f"Tip node index: {tip_node}, coordinates={coords[tip_node]}")
print(f"Computed mid deflection (y): {mid_disp_y:.6f}")
print(f"Analytical Euler-Bernoulli deflection: {w_analytical:.6f}")

# --- Evaluate error ---
error = abs(mid_disp_y - w_analytical)
print(f"Absolute error = {error:.6e}")

# --- Plot the mesh with the final deformed shape ---

img_name = "implemented_example.gif"
fname = str(tutorials_dir / img_name)
viz.make_deformation_gif(displacements_all, coords, connect, ele_type, fname)


