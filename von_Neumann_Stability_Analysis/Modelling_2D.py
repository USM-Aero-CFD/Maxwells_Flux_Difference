import math
import numpy as np
from numpy import linalg as LA

'''
    H_x = \hat{H}_x exp(- beta * \Delta x)
    H_y = \hat{H}_y exp(- beta * \Delta x)
    E_z = \hat{E}_z exp(- beta * \Delta x)

    | H_x |   | 0       0        0  |       |     0      0  1/mu |
    | H_y | = | 0       0     -1/mu | n_x + |     0      0   0   | n_y
    | E_z |   | 0  -1/epsilon    0  |       | 1/epsilon  0   0   |

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    | H_xj^n |                         | H_xi^n |
    | H_yj^n | = exp(-beta*\Delta x_j) | H_yi^n |
    | E_zj^n |                         | E_zi^n |

    | \Phi(H_x^T,n) |       | 0       0              0                                       | | H_xj^T,n |       |     0          0       (1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3mu | | H_xj^T,n |
    | \Phi(H_y^T,n) | = n_x | 0       0     -(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3mu | | H_yj^T,n | + n_y |     0          0                             0                         | | H_yj^T,n |
    | \Phi(E_z^T,n) |       | 0  -(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3epsilon    0  | | E_zj^T,n |       | (1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3epsilon        0      0   | | E_zj^T,n |

                        |                               0                                               0                                       n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3mu  | | H_xj^n |
                      = |                               0                                               0                                       -n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3mu | | H_yj^n |
                        | n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3epsilon   -n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3epsilon                            0                           | | E_zj^n |

    | k_1(H_x^T,n) |                 | \Phi(H_x^T,n) |
    | k_1(H_y^T,n) | = \sum_{T\cupi} | \Phi(H_y^T,n) |
    | k_1(E_z^T,n) |                 | \Phi(E_z^T,n) |

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    | H_xj^n+1/2 |   | H_xj^n |                                  |                                 0                                                         0                                   n_y(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/3mu   | | H_xj'^n |
    | H_yj^n+1/2 | = | H_yj^n | - \Delta t / 2 S_j \sum_{T\cupj} |                                 0                                                         0                                   -n_x(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/3mu  | | H_yj'^n |
    | E_zj^n+1/2 |   | E_zj^n |                                  | n_y(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/3epsilon   -n_x(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/3epsilon                               0                              | | E_zj'^n |

    | \Phi(H_x^T,n+1/2) |   |                               0                                               0                                       n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3mu  | | H_xj^n+1/2 |
    | \Phi(H_y^T,n+1/2) | = |                               0                                               0                                       -n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3mu | | H_yj^n+1/2 |
    | \Phi(E_z^T,n+1/2) |   | n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3epsilon   -n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3epsilon                            0                           | | E_zj^n+1/2 |

    | k_2(H_x^T,n+1/2) |                 | \Phi(H_x^T,n+1/2) |
    | k_2(H_y^T,n+1/2) | = \sum_{T\cupi} | \Phi(H_y^T,n+1/2) |
    | k_2(E_z^T,n+1/2) |                 | \Phi(E_z^T,n+1/2) |

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    | H_xj^n+1 |   | H_xj^n |                                |                               0                                                           0                                    n_y(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/3mu  | | H_xj'^n+1/2 |
    | H_yj^n+1 | = | H_yj^n | - \Delta t / S_j \sum_{T\cupj} |                               0                                                           0                                    -n_x(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/3mu | | H_yj'^n+1/2 |
    | E_zj^n+1 |   | E_zj^n |                                | n_y(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/3epsilon   -n_x(1+\sum_j'\inT,j'\notj exp(-beta*\Delta x_j'))/3epsilon                                 0                            | | E_zj'^n+1/2 |

    | \Phi(H_x^T,n+1) |   |                               0                                               0                                       n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3mu  | | H_xj^n+1 |
    | \Phi(H_y^T,n+1) | = |                               0                                               0                                       -n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3mu | | H_yj^n+1 |
    | \Phi(E_z^T,n+1) |   | n_y(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3epsilon   -n_x(1+\sum_j\inT,j\noti exp(-beta*\Delta x_j))/3epsilon                            0                           | | E_zj^n+1 |

    | k_3(H_x^T,n+1) |                 | \Phi(H_x^T,n+1) |
    | k_3(H_y^T,n+1) | = \sum_{T\cupi} | \Phi(H_y^T,n+1) |
    | k_3(E_z^T,n+1) |                 | \Phi(E_z^T,n+1) |

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    | H_xi^n+1 |   | H_xi^n |                   | k_1(H_x^T,n) + 4 k_2(H_x^T,n+1/2) + k_3(H_x^T,n+1) |
    | H_yi^n+1 | = | H_yi^n | - \Delta t / 6S_i | k_1(H_y^T,n) + 4 k_2(H_y^T,n+1/2) + k_3(H_y^T,n+1) |
    | E_zi^n+1 |   | E_zi^n |                   | k_1(E_z^T,n) + 4 k_2(E_z^T,n+1/2) + k_3(E_z^T,n+1) |
'''

# constant variables
Hx_0 = 1.0
Hy_0 = 1.0
Ez_0 = 1.0
beta = math.pi
omega = math.pi
mu = 1.0
epsilon = 1.0

def area_triangle(point_1, point_2, point_3):
    x1 = point_1['x']
    y1 = point_1['y']
    x2 = point_2['x']
    y2 = point_2['y']
    x3 = point_3['x']
    y3 = point_3['y']
    return abs((1.0 / 2.0) * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)))

def u(u_0, x, y, t):
    return complex(u_0 * math.cos(omega * t - beta * x - beta * y), u_0 * math.sin(omega * t - beta * x - beta * y))

def u_delta(u_0, x_delta, y_delta, t_delta):
    return complex(u_0 * math.cos(omega * t_delta - beta * x_delta - beta * y_delta), u_0 * math.sin(omega * t_delta - beta * x_delta - beta * y_delta))

def median(triangle):
    coordinate_1 = triangle[1]
    coordinate_2 = triangle[2]
    return {'nx': (coordinate_2['y'] - coordinate_1['y']) / 2.0, 'ny': - (coordinate_2['x'] - coordinate_1['x']) / 2.0} 

def u_average(triangle, u_0, t):
    coordinate_1 = triangle[0]
    coordinate_2 = triangle[1]
    coordinate_3 = triangle[2]
    u_1 = u(u_0, coordinate_1['x'], coordinate_1['y'], t)
    u_2 = u(u_0, coordinate_2['x'], coordinate_2['y'], t)
    u_3 = u(u_0, coordinate_3['x'], coordinate_3['y'], t)
    return (u_1 + u_2 + u_3) / 3.0

def flux_local_triangle(triangle, time_delta):
    normals = median(triangle)
    m_0 = u_average(triangle, Hx_0, time_delta)
    m_1 = u_average(triangle, Hy_0, time_delta)
    m_2 = u_average(triangle, Ez_0, time_delta)
    flux_matrix = np.matrix([[0 * m_0, 0 * m_1, m_2 * normals['ny'] / mu], [0 * m_0, 0 * m_1, - m_2 * normals['nx'] / mu], [m_0 * normals['ny'] / epsilon, - m_1 * normals['nx'] / epsilon, 0 * m_2]])
    
    return flux_matrix

def flux_nodal(triangles, time):
    sum = np.matrix([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

    for idx, triangle in enumerate(triangles):
        flux = flux_local_triangle(triangle, time)
        sum = np.add(sum, flux)
    
    return sum

'''
    G^n+1 U^n+1 = G^n U^n
    The eigenvalues \lambda could be obtained from
    det(G^n - \lambda G^n+1) = 0
'''
def nodal_update(time_delta, median_area, file_stream):
    #######################################################################################################
    # time level n
    #######################################################################################################
    # calculate flux at time level n
    time = 0.0
    k1 = flux_nodal(triangles, time)
    
    #######################################################################################################
    # time level n + 1/2
    #######################################################################################################
    # calculate flux at time level n + 1/2
    time = time_delta / 2.0
    k2 = flux_nodal(triangles, time)

    #######################################################################################################
    # Simpson's flux integration
    #######################################################################################################
    identity = np.matrix([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    ############################################  G^n  ########################################
    flux = np.add(k1, 4.0 * k2)
    flux = (time_delta / (6.0 * median_area)) * flux
    G_n = np.subtract(identity, flux)
    ###########################################  G^n+1  #######################################
    # time level n + 1
    ########################
    # calculate spatial flux only, will move the RHS flux terms to LHS: U^n+1
    time = 0.0
    k3 = flux_nodal(triangles, time)
    flux = (time_delta / (6.0 * median_area)) * k3
    G_n_1 = np.add(identity, flux)


    # print(f'Inverse of G^n+1 is {LA.inv(G_n_1)}')
    # print(f'Eigenvalues of (G^n+1)^-1 * G^n is {LA.eigvals(np.matmul(LA.inv(G_n_1), G_n))}')

    # inverse the G^n+1 matrix ---> get the amplification matrix ---> find its eigenvalues ---> get only the real part of eigenvalues
    eigenvalues_real = LA.eigvals(np.matmul(LA.inv(G_n_1), G_n)).real
    eigenvalues_abs = list(map(abs, eigenvalues_real))
    print(f'The largest eigenvalue is {max(eigenvalues_abs)}')
    file_stream.write("{:<30}".format(max(eigenvalues_abs)))
    file_stream.write("\n")

amplification_factor_2D = open('AmplificationFactor_2D_Analytic.txt', 'w')
for x_delta in np.arange(0.01, 1.0, 0.01):
    # hexagon stencil
    y_delta = x_delta
    vertex_0 = {'x': 0.0 * x_delta, 'y': 0.0 * y_delta}
    vertex_1 = {'x': 1.0 * x_delta, 'y': 0.0 * y_delta}
    vertex_2 = {'x': 0.5 * x_delta, 'y': (math.sqrt(3) / 2.0) * y_delta}
    vertex_3 = {'x': - 0.5 * x_delta, 'y': (math.sqrt(3) / 2.0) * y_delta}
    vertex_4 = {'x': - 1.0 * x_delta, 'y': 0.0 * y_delta}
    vertex_5 = {'x': - 0.5 * x_delta, 'y': (- math.sqrt(3) / 2.0) * y_delta}
    vertex_6 = {'x': 0.5 * x_delta, 'y': (- math.sqrt(3) / 2.0) * y_delta}
    vertices = [vertex_0, vertex_1, vertex_2, vertex_3, vertex_4, vertex_5, vertex_6]

    triangle_A = [vertex_0, vertex_1, vertex_2]
    triangle_B = [vertex_0, vertex_2, vertex_3]
    triangle_C = [vertex_0, vertex_3, vertex_4]
    triangle_D = [vertex_0, vertex_4, vertex_5]
    triangle_E = [vertex_0, vertex_5, vertex_6]
    triangle_F = [vertex_0, vertex_6, vertex_1]
    triangles = [triangle_A, triangle_B, triangle_C, triangle_D, triangle_E, triangle_F]

    # median cell area
    median_cell_area = 0.0
    for triangle in triangles:
        median_cell_area = median_cell_area + area_triangle(triangle[0], triangle[1], triangle[2])
    median_cell_area = median_cell_area / 3.0
    # print(f'Median dual cell area is: {median_cell_area}')

    # nodal update with fixed time_delta
    time_delta = 0.001
    print(f'x_delta is {x_delta}.')
    amplification_factor_2D.write("{:<10}".format(round(x_delta, 3)))
    nodal_update(time_delta, median_cell_area, amplification_factor_2D)
amplification_factor_2D.close()