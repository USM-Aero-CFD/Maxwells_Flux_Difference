import math
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg as LA

'''
    E_x = \hat{E}_x exp(- beta * \Delta x)
    E_y = \hat{E}_y exp(- beta * \Delta x)
    E_z = \hat{E}_z exp(- beta * \Delta x)
    H_x = \hat{H}_x exp(- beta * \Delta x) 
    H_y = \hat{H}_y exp(- beta * \Delta x)
    H_z = \hat{H}_z exp(- beta * \Delta x)

    | E_x |   | 0    0      0    0      0         0       |       |  0     0    0      0      0   -1/epsilon |       |  0      0    0        0      1/epsilon   0 |
    | E_y |   | 0    0      0    0      0       1/epsilon |       |  0     0    0      0      0       0      |       |  0      0    0    -1/epsilon     0       0 |
    | E_z | = | 0    0      0    0  -1/epsilon    0       | n_x + |  0     0    0   1/epsilon 0       0      | n_y + |  0      0    0        0          0       0 | n_z
    | H_x |   | 0    0      0    0      0         0       |       |  0     0   1/mu    0      0       0      |       |  0    -1/mu  0        0          0       0 |
    | H_y |   | 0    0    -1/mu  0      0         0       |       |  0     0    0      0      0       0      |       | 1/mu    0    0        0          0       0 |
    | H_z |   | 0   1/mu    0    0      0         0       |       | -1/mu  0    0      0      0       0      |       |  0      0    0        0          0       0 |

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    | E_xj^n |                         | E_xi^n |
    | E_yj^n |                         | E_yi^n |
    | E_zj^n | = exp(-beta*\Delta x_j) | E_zi^n |
    | H_xj^n |                         | H_xi^n |
    | H_yj^n |                         | H_yi^n |
    | H_zj^n |                         | H_zi^n |

    | \Phi(E_x^T,n) |   |                   0                                                           0                                                   0                                                       0                                       n_z(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon     -n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon    | | E_xj^n |
    | \Phi(E_y^T,n) |   |                   0                                                           0                                                   0                                  -n_z(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon                                 0                               n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon     | | E_yj^n |
    | \Phi(E_z^T,n) | = |                   0                                                           0                                                   0                                   n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon     -n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon                                0                               | | E_zj^n |
    | \Phi(H_x^T,n) |   |                   0                                   -n_z(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu     n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu                          0                                                                   0                                                           0                               | | H_xj^n |
    | \Phi(H_y^T,n) |   | n_z(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu                            0                               -n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu                         0                                                                   0                                                           0                               | | H_yj^n |
    | \Phi(H_z^T,n) |   | -n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu   n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu                          0                                                       0                                                                   0                                                           0                               | | H_zj^n |

    | k_1(E_x^T,n) |                 | \Phi(E_x^T,n) |
    | k_1(E_y^T,n) |                 | \Phi(E_y^T,n) |
    | k_1(E_z^T,n) | = \sum_{T\cupi} | \Phi(E_z^T,n) |    
    | k_1(H_x^T,n) |                 | \Phi(H_x^T,n) |
    | k_1(H_y^T,n) |                 | \Phi(H_y^T,n) |
    | k_1(H_z^T,n) |                 | \Phi(H_z^T,n) |

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    | E_xj^n+1/2 |   | E_xj^n |                                  |                   0                                                              0                                                           0                                                               0                                  n_z(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4epsilon     -n_y(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4epsilon   | | E_xj'^n |
    | E_yj^n+1/2 |   | E_yj^n |                                  |                   0                                                              0                                                           0                                  -n_z(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4epsilon                                      0                             n_x(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4epsilon    | | E_yj'^n |
    | E_zj^n+1/2 | = | E_zj^n | - \Delta t / 2 V_j \sum_{T\cupj} |                   0                                                              0                                                           0                                   n_y(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4epsilon     -n_x(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4epsilon                                  0                               | | E_zj'^n |
    | H_xj^n+1/2 |   | H_xj^n |                                  |                   0                                      -n_z(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4mu     n_y(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4mu                                    0                                                                   0                                                           0                               | | H_xj'^n |
    | H_yj^n+1/2 |   | H_yj^n |                                  | n_z(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4mu                            0                                   -n_x(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4mu                                  0                                                                   0                                                           0                               | | H_yj'^n |
    | H_zj^n+1/2 |   | H_zj^n |                                  | -n_y(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4mu   n_x(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4mu                               0                                                               0                                                                   0                                                           0                               | | H_zj'^n |


    | \Phi(E_x^T,n+1/2) |   |                   0                                                           0                                                   0                                                       0                                       n_z(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon     -n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon    | | E_xj^n+1/2 |
    | \Phi(E_y^T,n+1/2) |   |                   0                                                           0                                                   0                                  -n_z(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon                                 0                               n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon     | | E_yj^n+1/2 |
    | \Phi(E_z^T,n+1/2) | = |                   0                                                           0                                                   0                                   n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon     -n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon                                0                               | | E_zj^n+1/2 |
    | \Phi(H_x^T,n+1/2) |   |                   0                                   -n_z(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu     n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu                          0                                                                   0                                                           0                               | | H_xj^n+1/2 |
    | \Phi(H_y^T,n+1/2) |   | n_z(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu                            0                               -n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu                         0                                                                   0                                                           0                               | | H_yj^n+1/2 |
    | \Phi(H_z^T,n+1/2) |   | -n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu   n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu                          0                                                       0                                                                   0                                                           0                               | | H_zj^n+1/2 |

    | k_2(E_x^T,n+1/2) |                 | \Phi(E_x^T,n+1/2) |
    | k_2(E_y^T,n+1/2) |                 | \Phi(E_y^T,n+1/2) |
    | k_2(E_z^T,n+1/2) | = \sum_{T\cupi} | \Phi(E_z^T,n+1/2) |    
    | k_2(H_x^T,n+1/2) |                 | \Phi(H_x^T,n+1/2) |
    | k_2(H_y^T,n+1/2) |                 | \Phi(H_y^T,n+1/2) |
    | k_2(H_z^T,n+1/2) |                 | \Phi(H_z^T,n+1/2) |

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    | E_xj^n+1 |   | E_xj^n |                                |                   0                                                              0                                                           0                                                               0                                  n_z(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4epsilon     -n_y(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4epsilon   | | E_xj'^n+1/2 |
    | E_yj^n+1 |   | E_yj^n |                                |                   0                                                              0                                                           0                                  -n_z(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4epsilon                                      0                             n_x(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4epsilon    | | E_yj'^n+1/2 |
    | E_zj^n+1 | = | E_zj^n | - \Delta t / V_j \sum_{T\cupj} |                   0                                                              0                                                           0                                   n_y(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4epsilon     -n_x(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4epsilon                                  0                               | | E_zj'^n+1/2 |
    | H_xj^n+1 |   | H_xj^n |                                |                   0                                      -n_z(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4mu     n_y(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4mu                                    0                                                                   0                                                           0                               | | H_xj'^n+1/2 |
    | H_yj^n+1 |   | H_yj^n |                                | n_z(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4mu                            0                                   -n_x(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4mu                                  0                                                                   0                                                           0                               | | H_yj'^n+1/2 |
    | H_zj^n+1 |   | H_zj^n |                                | -n_y(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4mu   n_x(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/4mu                               0                                                               0                                                                   0                                                           0                               | | H_zj'^n+1/2 |


    | \Phi(E_x^T,n+1) |   |                   0                                                           0                                                   0                                                       0                                       n_z(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon     -n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon    | | E_xj^n+1 |
    | \Phi(E_y^T,n+1) |   |                   0                                                           0                                                   0                                  -n_z(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon                                 0                               n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon     | | E_yj^n+1 |
    | \Phi(E_z^T,n+1) | = |                   0                                                           0                                                   0                                   n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon     -n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4epsilon                                0                               | | E_zj^n+1 |
    | \Phi(H_x^T,n+1) |   |                   0                                   -n_z(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu     n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu                          0                                                                   0                                                           0                               | | H_xj^n+1 |
    | \Phi(H_y^T,n+1) |   | n_z(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu                            0                               -n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu                         0                                                                   0                                                           0                               | | H_yj^n+1 |
    | \Phi(H_z^T,n+1) |   | -n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu   n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/4mu                          0                                                       0                                                                   0                                                           0                               | | H_zj^n+1 |

    | k_2(E_x^T,n+1) |                 | \Phi(E_x^T,n+1) |
    | k_2(E_y^T,n+1) |                 | \Phi(E_y^T,n+1) |
    | k_2(E_z^T,n+1) | = \sum_{T\cupi} | \Phi(E_z^T,n+1) |    
    | k_2(H_x^T,n+1) |                 | \Phi(H_x^T,n+1) |
    | k_2(H_y^T,n+1) |                 | \Phi(H_y^T,n+1) |
    | k_2(H_z^T,n+1) |                 | \Phi(H_z^T,n+1) |

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    | E_xi^n+1 |   | E_xi^n |                   | k_1(E_x^T,n) + 4 k_2(E_x^T,n+1/2) + k_3(E_x^T,n+1) |
    | E_yi^n+1 |   | E_yi^n |                   | k_1(E_y^T,n) + 4 k_2(E_y^T,n+1/2) + k_3(E_y^T,n+1) |
    | E_zi^n+1 | = | E_zi^n | - \Delta t / 6V_i | k_1(E_z^T,n) + 4 k_2(E_z^T,n+1/2) + k_3(E_z^T,n+1) |
    | H_xi^n+1 |   | H_xi^n |                   | k_1(H_x^T,n) + 4 k_2(H_x^T,n+1/2) + k_3(H_x^T,n+1) |
    | H_yi^n+1 |   | H_yi^n |                   | k_1(H_y^T,n) + 4 k_2(H_y^T,n+1/2) + k_3(H_y^T,n+1) |
    | H_zi^n+1 |   | H_zi^n |                   | k_1(H_z^T,n) + 4 k_2(H_z^T,n+1/2) + k_3(H_z^T,n+1) |      
'''

# constant variables
Ex_0 = 1.0
Ey_0 = 1.0
Ez_0 = 1.0
Hx_0 = 1.0
Hy_0 = 1.0
Hz_0 = 1.0
beta = math.pi
omega = math.pi
mu = 1.0
epsilon = 1.0

def u(u_0, x, y, z, t):
    return complex(u_0 * math.cos(omega * t - beta * x - beta * y - beta * z), u_0 * math.sin(omega * t - beta * x - beta * y - beta * z))

def u_delta(u_0, x_delta, y_delta, z_delta, t_delta):
    return complex(u_0 * math.cos(omega * t_delta - beta * x_delta - beta * y_delta - beta * z_delta), u_0 * math.sin(omega * t_delta - beta * x_delta - beta * y_delta - beta * z_delta))

def displacement(vector_1, vector_2):
    x1 = vector_1['x']
    y1 = vector_1['y']
    z1 = vector_1['z']
    x2 = vector_2['x']
    y2 = vector_2['y']
    z2 = vector_2['z']
    return {'x': (x2 - x1), 'y': (y2 - y1), 'z': (z2 - z1)}

def dot_product(vector_1, vector_2):
    x1 = vector_1['x']
    y1 = vector_1['y']
    z1 = vector_1['z']
    x2 = vector_2['x']
    y2 = vector_2['y']
    z2 = vector_2['z']
    return x1 * x2 + y1 * y2 + z1 * z2

def cross_product(vector_1, vector_2):
    x1 = vector_1['x']
    y1 = vector_1['y']
    z1 = vector_1['z']
    x2 = vector_2['x']
    y2 = vector_2['y']
    z2 = vector_2['z']
    return {'x': (y1 * z2 - y2 * z1), 'y': (x2 * z1 - x1 * z2), 'z': (x1 * y2 - x2 * y1)}

def median(tetrahedron):
    coordinate_1 = tetrahedron[1]
    coordinate_2 = tetrahedron[2]
    coordinate_3 = tetrahedron[3]
    edge_1 = displacement(coordinate_2, coordinate_3)
    edge_2 = displacement(coordinate_2, coordinate_1)
    opposite_plane = cross_product(edge_1, edge_2)
    return {'nx': - opposite_plane['x'] / 3.0, 'ny': - opposite_plane['y'] / 3.0, 'nz': - opposite_plane['z'] / 3.0}

def volume_tetrahedron(point_1, point_2, point_3, point_4):
    edge_1 = displacement(point_3, point_4)
    edge_2 = displacement(point_3, point_2)
    opposite_plane = cross_product(edge_1, edge_2)
    height = displacement(point_3, point_1)
    return (1.0 / 3.0) * abs(dot_product(height, opposite_plane))

def u_average(tetrahedron, u_0, t):
    coordinate_1 = tetrahedron[0]
    coordinate_2 = tetrahedron[1]
    coordinate_3 = tetrahedron[2]
    coordinate_4 = tetrahedron[3]
    u_1 = u(u_0, coordinate_1['x'], coordinate_1['y'], coordinate_1['z'], t)
    u_2 = u(u_0, coordinate_2['x'], coordinate_2['y'], coordinate_2['z'], t)
    u_3 = u(u_0, coordinate_3['x'], coordinate_3['y'], coordinate_3['z'], t)
    u_4 = u(u_0, coordinate_4['x'], coordinate_4['y'], coordinate_4['z'], t)
    return (u_1 + u_2 + u_3 + u_4) / 4.0

def flux_local_tetrahedron(tetrahedron, time_delta):
    normals = median(tetrahedron)
    m_0 = u_average(tetrahedron, Ex_0, time_delta)
    m_1 = u_average(tetrahedron, Ey_0, time_delta)
    m_2 = u_average(tetrahedron, Ez_0, time_delta)
    m_3 = u_average(tetrahedron, Hx_0, time_delta)
    m_4 = u_average(tetrahedron, Hy_0, time_delta)
    m_5 = u_average(tetrahedron, Hz_0, time_delta)
    flux_matrix = np.matrix([[0 * m_0, 0 * m_1, 0 * m_2, 0 * m_3, m_4 * normals['nz'] / epsilon, - m_5 * normals['ny'] / epsilon], 
        [0 * m_0, 0 * m_1, 0 * m_2, - m_3 * normals['nz'] / epsilon, 0 * m_4, m_5 * normals['nx'] / epsilon], 
        [0 * m_0, 0 * m_1, 0 * m_2, m_3 * normals['ny'] / epsilon, - m_4 * normals['nx'] / epsilon, 0 * m_5], 
        [0 * m_0, - m_1 * normals['nz'] / mu, m_2 * normals['ny'] / mu, 0 * m_3, 0 * m_4, 0 * m_5], 
        [m_0 * normals['nz'] / mu, 0 * m_1, - m_2 * normals['nx'] / mu, 0 * m_3, 0 * m_4, 0 * m_5], 
        [- m_0 * normals['ny'] / mu, m_1 * normals['nx'] / mu, 0 * m_2, 0 * m_3, 0 * m_4, 0 * m_5]])

    return flux_matrix

def flux_nodal(tetrahedra, time):
    sum = np.matrix([[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
    
    for idx, tetrahedron in enumerate(tetrahedra):
        flux = flux_local_tetrahedron(tetrahedron, time)
        sum = np.add(sum, flux)
    
    return sum

'''
    G^n+1 U^n+1 = G^n U^n
    The eigenvalues \lambda could be obtained from
    det(G^n - \lambda G^n+1) = 0
'''
def nodal_update(time_delta, median_cell_volume, file_stream):
    #######################################################################################################
    # time level n
    #######################################################################################################
    # calculate flux at time level n
    time = 0.0
    k1 = flux_nodal(tetrahedra, time)

    #######################################################################################################
    # time level n + 1/2
    #######################################################################################################
    # calculate flux at time level n + 1/2
    time = time_delta / 2.0
    k2 = flux_nodal(tetrahedra, time)

    #######################################################################################################
    # Simpson's flux integration
    #######################################################################################################
    identity = np.matrix([[1.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
        [0.0, 1.0, 0.0, 0.0, 0.0, 0.0], 
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0], 
        [0.0, 0.0, 0.0, 1.0, 0.0, 0.0], 
        [0.0, 0.0, 0.0, 0.0, 1.0, 0.0], 
        [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]])
    ############################################  G^n  ########################################
    flux = np.add(k1, 4.0 * k2)
    flux = (time_delta / (6.0 * median_cell_volume)) * flux
    G_n = np.subtract(identity, flux)
    ###########################################  G^n+1  #######################################
    # time level n + 1
    ########################
    # calculate spatial flux only, will move the RHS flux terms to LHS: U^n+1
    time = 0.0
    k3 = flux_nodal(tetrahedra, time)
    flux = (time_delta / (6.0 * median_cell_volume)) * k3
    G_n_1 = np.add(identity, flux)

    # print(f'Inverse of G^n+1 is {LA.inv(G_n_1)}')
    # print(f'Eigenvalues of (G^n+1)^-1 * G^n is {LA.eigvals(np.matmul(LA.inv(G_n_1), G_n))}')

    # inverse the G^n+1 matrix ---> get the amplification matrix ---> find its eigenvalues ---> get only the real part of eigenvalues
    eigenvalues_real = LA.eigvals(np.matmul(LA.inv(G_n_1), G_n)).real
    eigenvalues_abs = list(map(abs, eigenvalues_real))
    print(f'The largest eigenvalue is {max(eigenvalues_abs)}')
    file_stream.write("{:<30}".format(max(eigenvalues_abs)))
    file_stream.write("\n")

amplification_factor_3D = open('AmplificationFactor_3D_Analytic.txt', 'w')
for x_delta in np.arange(0.01, 1.0, 0.01):
    # icosahedron stencil
    y_delta = z_delta = x_delta
    vertex_0 = {'x': 0.0 * x_delta, 'y': 0.0 * y_delta, 'z': 0.0 * z_delta}
    vertex_1 = {'x': 1.0 * x_delta, 'y': 0.0 * y_delta, 'z': 0.0 * z_delta}
    vertex_2 = {'x': - 1.0 * x_delta, 'y': 0.0 * y_delta, 'z': 0.0 * z_delta}
    vertex_3 = {'x': (1.0 / math.sqrt(5.0)) * x_delta, 'y': (2.0 / math.sqrt(5.0)) * y_delta, 'z': 0.0 * z_delta}
    vertex_4 = {'x': (- 1.0 / math.sqrt(5.0)) * x_delta, 'y': (- 2.0 / math.sqrt(5.0)) * y_delta, 'z': 0.0 * z_delta}
    vertex_5 = {'x': (1.0 / math.sqrt(5.0)) * x_delta, 'y': ((1.0 - 1.0 / math.sqrt(5.0)) / 2.0) * y_delta, 'z': (math.sqrt((1.0 + 1.0 / math.sqrt(5.0)) / 2.0)) * z_delta}
    vertex_6 = {'x': (1.0 / math.sqrt(5.0)) * x_delta, 'y': ((1.0 - 1.0 / math.sqrt(5.0)) / 2.0) * y_delta, 'z': (- math.sqrt((1.0 + 1.0 / math.sqrt(5.0)) / 2.0)) * z_delta}
    vertex_7 = {'x': (- 1.0 / math.sqrt(5.0)) * x_delta, 'y': (- (1.0 - 1.0 / math.sqrt(5.0)) / 2.0) * y_delta, 'z': (- math.sqrt((1.0 + 1.0 / math.sqrt(5.0)) / 2.0)) * z_delta}
    vertex_8 = {'x': (- 1.0 / math.sqrt(5.0)) * x_delta, 'y': (- (1.0 - 1.0 / math.sqrt(5.0)) / 2.0) * y_delta, 'z': (math.sqrt((1.0 + 1.0 / math.sqrt(5.0)) / 2.0)) * z_delta}
    vertex_9 = {'x': (1.0 / math.sqrt(5.0)) * x_delta, 'y': ((- 1.0 - 1.0 / math.sqrt(5.0)) / 2.0) * y_delta, 'z': (math.sqrt((1.0 - 1.0 / math.sqrt(5.0)) / 2.0)) * z_delta}
    vertex_10 = {'x': (1.0 / math.sqrt(5.0)) * x_delta, 'y': ((- 1.0 - 1.0 / math.sqrt(5.0)) / 2.0) * y_delta, 'z': (- math.sqrt((1.0 - 1.0 / math.sqrt(5.0)) / 2.0)) * z_delta}
    vertex_11 = {'x': (- 1.0 / math.sqrt(5.0)) * x_delta, 'y': (- (- 1.0 - 1.0 / math.sqrt(5.0)) / 2.0) * y_delta, 'z': (- math.sqrt((1.0 - 1.0 / math.sqrt(5.0)) / 2.0)) * z_delta}
    vertex_12 = {'x': (- 1.0 / math.sqrt(5.0)) * x_delta, 'y': (- (- 1.0 - 1.0 / math.sqrt(5.0)) / 2.0) * y_delta, 'z': (math.sqrt((1.0 - 1.0 / math.sqrt(5.0)) / 2.0)) * z_delta}
    vertices = [vertex_0, vertex_1, vertex_2, vertex_3, vertex_4, vertex_5, vertex_6, vertex_7, vertex_8, vertex_9, vertex_10, vertex_11, vertex_12]

    #  local vertices of tetrahedra
    tetrahedron_A = [vertex_0, vertex_1, vertex_5, vertex_3]
    tetrahedron_B = [vertex_0, vertex_1, vertex_3, vertex_6]
    tetrahedron_C = [vertex_0, vertex_1, vertex_6, vertex_10]
    tetrahedron_D = [vertex_0, vertex_1, vertex_9, vertex_5]
    tetrahedron_E = [vertex_0, vertex_1, vertex_10, vertex_9]
    tetrahedron_F = [vertex_0, vertex_2, vertex_4, vertex_7]
    tetrahedron_G = [vertex_0, vertex_2, vertex_7, vertex_11]
    tetrahedron_H = [vertex_0, vertex_2, vertex_8, vertex_4]
    tetrahedron_I = [vertex_0, vertex_2, vertex_12, vertex_8]
    tetrahedron_J = [vertex_0, vertex_2, vertex_11, vertex_12]
    tetrahedron_K = [vertex_0, vertex_3, vertex_5, vertex_12]
    tetrahedron_L = [vertex_0, vertex_3, vertex_11, vertex_6]
    tetrahedron_M = [vertex_0, vertex_3, vertex_12, vertex_11]
    tetrahedron_N = [vertex_0, vertex_4, vertex_8, vertex_9]
    tetrahedron_O = [vertex_0, vertex_4, vertex_9, vertex_10]
    tetrahedron_P = [vertex_0, vertex_4, vertex_10, vertex_7]
    tetrahedron_Q = [vertex_0, vertex_5, vertex_8, vertex_12]
    tetrahedron_R = [vertex_0, vertex_5, vertex_9, vertex_8]
    tetrahedron_S = [vertex_0, vertex_6, vertex_7, vertex_10]
    tetrahedron_T = [vertex_0, vertex_6, vertex_11, vertex_7]
    tetrahedra = [tetrahedron_A, tetrahedron_B, tetrahedron_C, tetrahedron_D, tetrahedron_E, 
                    tetrahedron_F, tetrahedron_G, tetrahedron_H, tetrahedron_I, tetrahedron_J,
                    tetrahedron_K, tetrahedron_L, tetrahedron_M, tetrahedron_N, tetrahedron_O,
                    tetrahedron_P, tetrahedron_Q, tetrahedron_R, tetrahedron_S, tetrahedron_T]

    # median cell area
    median_cell_volume = 0.0
    for tetrahedron in tetrahedra:
        median_cell_volume = median_cell_volume + volume_tetrahedron(tetrahedron[0], tetrahedron[1], tetrahedron[2], tetrahedron[3])
    median_cell_volume = median_cell_volume / 4.0
    # print(f'Median dual cell volume is: {median_cell_volume}')

    # nodal update with fixed time_delta
    time_delta = 0.001
    print(f'x_delta is {x_delta}.')
    amplification_factor_3D.write("{:<10}".format(round(x_delta, 3)))
    nodal_update(time_delta, median_cell_volume, amplification_factor_3D)
amplification_factor_3D.close()

############################################################################################
## visualizing all the vertices
############################################################################################
# # icosahedron stencil
# vertex_0 = {'x': 0.0, 'y': 0.0, 'z': 0.0}
# vertex_1 = {'x': 1.0, 'y': 0.0, 'z': 0.0}
# vertex_2 = {'x': - 1.0, 'y': 0.0, 'z': 0.0}
# vertex_3 = {'x': (1.0 / math.sqrt(5.0)), 'y': (2.0 / math.sqrt(5.0)), 'z': 0.0}
# vertex_4 = {'x': (- 1.0 / math.sqrt(5.0)), 'y': (- 2.0 / math.sqrt(5.0)), 'z': 0.0}
# vertex_5 = {'x': (1.0 / math.sqrt(5.0)), 'y': ((1.0 - 1.0 / math.sqrt(5.0)) / 2.0), 'z': (math.sqrt((1.0 + 1.0 / math.sqrt(5.0)) / 2.0))}
# vertex_6 = {'x': (1.0 / math.sqrt(5.0)), 'y': ((1.0 - 1.0 / math.sqrt(5.0)) / 2.0), 'z': (- math.sqrt((1.0 + 1.0 / math.sqrt(5.0)) / 2.0))}
# vertex_7 = {'x': (- 1.0 / math.sqrt(5.0)), 'y': (- (1.0 - 1.0 / math.sqrt(5.0)) / 2.0), 'z': (- math.sqrt((1.0 + 1.0 / math.sqrt(5.0)) / 2.0))}
# vertex_8 = {'x': (- 1.0 / math.sqrt(5.0)), 'y': (- (1.0 - 1.0 / math.sqrt(5.0)) / 2.0), 'z': (math.sqrt((1.0 + 1.0 / math.sqrt(5.0)) / 2.0))}
# vertex_9 = {'x': (1.0 / math.sqrt(5.0)), 'y': ((- 1.0 - 1.0 / math.sqrt(5.0)) / 2.0), 'z': (math.sqrt((1.0 - 1.0 / math.sqrt(5.0)) / 2.0))}
# vertex_10 = {'x': (1.0 / math.sqrt(5.0)), 'y': ((- 1.0 - 1.0 / math.sqrt(5.0)) / 2.0), 'z': (- math.sqrt((1.0 - 1.0 / math.sqrt(5.0)) / 2.0))}
# vertex_11 = {'x': (- 1.0 / math.sqrt(5.0)), 'y': (- (- 1.0 - 1.0 / math.sqrt(5.0)) / 2.0), 'z': (- math.sqrt((1.0 - 1.0 / math.sqrt(5.0)) / 2.0))}
# vertex_12 = {'x': (- 1.0 / math.sqrt(5.0)), 'y': (- (- 1.0 - 1.0 / math.sqrt(5.0)) / 2.0), 'z': (math.sqrt((1.0 - 1.0 / math.sqrt(5.0)) / 2.0))}
# vertices = [vertex_0, vertex_1, vertex_2, vertex_3, vertex_4, vertex_5, vertex_6, vertex_7, vertex_8, vertex_9, vertex_10, vertex_11, vertex_12]

# ax = plt.figure().add_subplot(projection='3d')
# for idx, vertex in enumerate(vertices):
#     x = vertex['x']
#     y = vertex['y']
#     z = vertex['z']
#     print('(%d, %f, %f, %f)'%(idx, x, y, z))
#     label = '%d' % (idx)
#     ax.text(x, y, z, label)

# ax.set_xlim(-1, 1)
# ax.set_ylim(-1, 1)
# ax.set_zlim(-1, 1)
# ax.set_xlabel('X axis')
# ax.set_ylabel('Y axis')
# ax.set_zlabel('Z axis')

# plt.show()


############################################################################################
## visualizing all the surfaces with their outward normals
############################################################################################
# median_0_A = median(vertex_1, vertex_5, vertex_3)
# median_0_B = median(vertex_1, vertex_3, vertex_6)
# median_0_C = median(vertex_1, vertex_6, vertex_10)
# median_0_D = median(vertex_1, vertex_9, vertex_5)
# median_0_E = median(vertex_1, vertex_10, vertex_9)
# median_0_F = median(vertex_2, vertex_4, vertex_7)
# median_0_G = median(vertex_2, vertex_7, vertex_11)
# median_0_H = median(vertex_2, vertex_8, vertex_4)
# median_0_I = median(vertex_2, vertex_12, vertex_8)
# median_0_J = median(vertex_2, vertex_11, vertex_12)
# median_0_K = median(vertex_3, vertex_5, vertex_12)
# median_0_L = median(vertex_3, vertex_11, vertex_6)
# median_0_M = median(vertex_3, vertex_12, vertex_11)
# median_0_N = median(vertex_4, vertex_8, vertex_9)
# median_0_O = median(vertex_4, vertex_9, vertex_10)
# median_0_P = median(vertex_4, vertex_10, vertex_7)
# median_0_Q = median(vertex_5, vertex_8, vertex_12)
# median_0_R = median(vertex_5, vertex_9, vertex_8)
# median_0_S = median(vertex_6, vertex_7, vertex_10)
# median_0_T = median(vertex_6, vertex_11, vertex_7)
# median_dual_normals = [median_0_A, median_0_B, median_0_C, median_0_D, median_0_E, median_0_F, median_0_G, median_0_H, median_0_I, median_0_J, 
#                         median_0_K, median_0_L, median_0_M, median_0_N, median_0_O, median_0_P, median_0_Q, median_0_R, median_0_S, median_0_T]

# def centroid(point_1, point_2, point_3):
#     return {'x': (point_1['x'] + point_2['x'] + point_3['x']) / 3.0, 
#             'y': (point_1['y'] + point_2['y'] + point_3['y']) / 3.0, 
#             'z': (point_1['z'] + point_2['z'] + point_3['z']) / 3.0 }

# plane_centroid_A = centroid(vertex_1, vertex_5, vertex_3)
# plane_centroid_B = centroid(vertex_1, vertex_3, vertex_6)
# plane_centroid_C = centroid(vertex_1, vertex_6, vertex_10)
# plane_centroid_D = centroid(vertex_1, vertex_9, vertex_5)
# plane_centroid_E = centroid(vertex_1, vertex_10, vertex_9)
# plane_centroid_F = centroid(vertex_2, vertex_4, vertex_7)
# plane_centroid_G = centroid(vertex_2, vertex_7, vertex_11)
# plane_centroid_H = centroid(vertex_2, vertex_8, vertex_4)
# plane_centroid_I = centroid(vertex_2, vertex_12, vertex_8)
# plane_centroid_J = centroid(vertex_2, vertex_11, vertex_12)
# plane_centroid_K = centroid(vertex_3, vertex_5, vertex_12)
# plane_centroid_L = centroid(vertex_3, vertex_11, vertex_6)
# plane_centroid_M = centroid(vertex_3, vertex_12, vertex_11)
# plane_centroid_N = centroid(vertex_4, vertex_8, vertex_9)
# plane_centroid_O = centroid(vertex_4, vertex_9, vertex_10)
# plane_centroid_P = centroid(vertex_4, vertex_10, vertex_7)
# plane_centroid_Q = centroid(vertex_5, vertex_8, vertex_12)
# plane_centroid_R = centroid(vertex_5, vertex_9, vertex_8)
# plane_centroid_S = centroid(vertex_6, vertex_7, vertex_10)
# plane_centroid_T = centroid(vertex_6, vertex_11, vertex_7)
# plane_centroids = [plane_centroid_A, plane_centroid_B, plane_centroid_C, plane_centroid_D, plane_centroid_E, plane_centroid_F, plane_centroid_G, plane_centroid_H, plane_centroid_I, plane_centroid_J, 
#                         plane_centroid_K, plane_centroid_L, plane_centroid_M, plane_centroid_N, plane_centroid_O, plane_centroid_P, plane_centroid_Q, plane_centroid_R, plane_centroid_S, plane_centroid_T]

# ax = plt.figure().add_subplot(projection='3d')

# # Make the grid
# x, y, z = [], [], []
# for centroid in plane_centroids:
#     x.append(centroid['x'])
#     y.append(centroid['y'])
#     z.append(centroid['z'])

# # Make the direction data for the arrows
# u, v, w = [], [], []
# for median_dual in median_dual_normals:
#     u.append(median_dual['nx'])
#     v.append(median_dual['ny'])
#     w.append(median_dual['nz'])

# ax.quiver(x, y, z, u, v, w, length=0.3, normalize=True)

# plt.show()