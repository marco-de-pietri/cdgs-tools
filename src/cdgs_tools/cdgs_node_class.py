"""
Module for the CDGS node class.
"""
import math
import pyvista as pv
from vtk import vtkSelectEnclosedPoints
import numpy as np

def write_vec(vec, npos, digits, spaces=1, spaces_new_line=0):
    """
    Line formatting with a specified number of spaces between numbers.
    """
    j = -1
    new_line_empty = '' * spaces_new_line
    line = new_line_empty
    space_str = ' ' * spaces
    for el_val in vec:
        j = j + 1
        if j == npos:
            line = line + '\n' + new_line_empty
            j = 0
        line = line + f"{el_val:.{digits}e}{space_str}"
    return line

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a_1 = math.cos(theta / 2.0)
    b_1, c_1, d_1 = -axis * math.sin(theta / 2.0)
    a_a, b_b, c_c, d_d = a_1 * a_1, b_1 * b_1, c_1 * c_1, d_1 * d_1
    b_c, a_d, a_c, a_b, b_d, c_d = ( b_1 * c_1,
                                     a_1 * d_1,
                                     a_1 * c_1,
                                     a_1 * b_1,
                                     b_1 * d_1,
                                     c_1 * d_1)
    res = np.array([[a_a + b_b - c_c - d_d, 2 * (b_c + a_d), 2 * (b_d - a_c)],
                     [2 * (b_c - a_d), a_a + c_c - b_b - d_d, 2 * (c_d + a_b)],
                     [2 * (b_d + a_c), 2 * (c_d - a_b), a_a + d_d - b_b - c_c]])
    return res

def extend_pnt(point,axis, distance):
    """
    this function returns a point by extending a point, along a vector,
    by a given distance
    """

    newpoint = [point[0] + axis[0]*distance,
                point[1] + axis[1]*distance,
                point[2] + axis[2]*distance]

    return newpoint


class CDGSNode:
    """
    CDGS node generic class.
    """
    def __init__(self) -> None:
        self.mesh_id = 0
        self.energy_bins = []
        self.energy_bin_vals = []
        self.energy_bin_errs = []
        self.node_intensity_s = 0
        self.node_intensity_m3_s = 0
        self.node_comment = ""
        self.cooling_time = 0
        self.cdgs_text = ""

        self.vtk_points = []
        self.vtk_points_text = ""

        self.vtk_cell_type = 0
        self.vtk_point_ids = []
        self.vtk_point_ids_text = ""


    def __str__(self) -> str:
        return (
            f"CDGSNode object with mesh id {self.mesh_id} "
            f"and node intensity {self.node_intensity_s:.2e}"
        )
class CDGSCylinderNode(CDGSNode):
    """
    class for cdgs cylinder node
    """

    def __init__(self) -> None:
        super().__init__()
        self.basis_cm = []
        self.end_point_cm = []
        self.axis_1 = []
        self.axis_2 = []
        self.radius_cm = 0
        self.height_cm = 0
        self.node_volume_cm3 = 0
        self.node_volume_m3 = 0
        self.lateral_surface_area_cm2 = 0



    def calculate_end_point(self) -> None:
        """
        Calculates the end point of the cylinder.
        """
        self.end_point_cm = [
            self.basis_cm[0] + self.axis_1[0] * self.height_cm,
            self.basis_cm[1] + self.axis_1[1] * self.height_cm,
            self.basis_cm[2] + self.axis_1[2] * self.height_cm,
        ]

    def generate_cdgs_text(self) -> None:
        """
        Generates the CDGS text.
        """
        text = ""
        text += f"mesh_id {self.mesh_id}\n"
        text += f"{self.node_comment}\n"
        text += f"Cooling_time {self.cooling_time:.6e}\n"
        text += f"total_source {self.node_intensity_s:.9e}\n"
        text += f"energy_type bins\n"
        text += f"energy_boundaries {len(self.energy_bins):d}\n"
        text += write_vec(self.energy_bins, 6, 7)
        text += "\n"
        text += f"mesh_type cyl\n"
        text += f"mesh_boundaries  2 2 2\n"
        text += write_vec(self.basis_cm, 6, 6, 2)
        text += "\n"
        text += write_vec(self.axis_1, 6, 6, 2)
        text += "\n"
        text += write_vec(self.axis_2, 6, 6, 2)
        text += "\n"
        text += f"      0.0000  {self.radius_cm:>11.4f}\n"
        text += f"      0.0000       1.0000\n"
        text += f"      0.0000  {self.height_cm:>11.4f}\n"
        text += f"source_data\n"
        text += f"1  {self.node_intensity_s:.8e}  {self.node_volume_cm3:.8e} 1\n"
        text += f"{self.mesh_id:d} 1.00000000e+00  {self.node_intensity_s:.8e}\n"
        text += write_vec(self.energy_bin_vals, 6, 4, 3, 2)
        text += "\n"
        text += write_vec(self.energy_bin_errs, 6, 4, 3, 2)
        text += "\n"

        text += f"end_source_data\n"
        self.cdgs_text = text

        return


    def calculate_volume(self) -> None:
        """
        Calculates the volume of the cylinder.
        """
        self.node_volume_cm3 = math.pi * self.height_cm * (self.radius_cm ** 2)
        self.node_volume_m3 =  self.node_volume_cm3 / 1e6


    def calculate_lateral_surface_area(self) -> None:
        """
        Calculates the lateral surface area of the cylinder.
        """
        self.lateral_surface_area_cm2 = 2 * math.pi * self.radius_cm * self.height_cm

    def generate_vtk_points(self) -> None:
        """
        Calculates the VTK points of the cylinder.
        """
        pipe_axis_1 = [self.axis_1[0], self.axis_1[1], self.axis_1[2]]
        pipe_axis_2 = [self.axis_2[0], self.axis_2[1], self.axis_2[2]]
        radius_m = self.radius_cm/100
        length_m = self.height_cm/100
        pipe_origin_m = [self.basis_cm[0]/100,
                         self.basis_cm[1]/100,
                         self.basis_cm[2]/100]

        pipe_end_m = extend_pnt(pipe_origin_m,
                                       pipe_axis_1,
                                       length_m)
        pipe_mid_m = extend_pnt(pipe_origin_m,
                                       pipe_axis_1,
                                       length_m/2)

        vector0   = pipe_axis_2

        vector45  = np.dot(rotation_matrix(pipe_axis_1,math.pi/4), vector0)
        vector90  = np.dot(rotation_matrix(pipe_axis_1,math.pi/2), vector0)
        vector135 = np.dot(rotation_matrix(pipe_axis_1,math.pi*3/4),vector0)
        vector180 = np.dot(rotation_matrix(pipe_axis_1,math.pi), vector0)
        vector225 = np.dot(rotation_matrix(pipe_axis_1,math.pi*5/4), vector0)
        vector270 = np.dot(rotation_matrix(pipe_axis_1,math.pi + math.pi/2),
                           vector0)
        vector315 = np.dot(rotation_matrix(pipe_axis_1,math.pi*7/4), vector0)

        # reference :   https://vtk.org/doc/release/5.2/html/a00103.html

        point00 = extend_pnt(pipe_origin_m,vector0,  radius_m)
        point01 = extend_pnt(pipe_origin_m,vector90, radius_m)
        point02 = extend_pnt(pipe_origin_m,vector180,radius_m)
        point03 = extend_pnt(pipe_origin_m,vector270,radius_m)
        point04 = extend_pnt(pipe_end_m,vector0,  radius_m)
        point05 = extend_pnt(pipe_end_m,vector90, radius_m)
        point06 = extend_pnt(pipe_end_m,vector180,radius_m)
        point07 = extend_pnt(pipe_end_m,vector270,radius_m)
        point08 = extend_pnt(pipe_origin_m,vector45,radius_m)
        point09 = extend_pnt(pipe_origin_m,vector135,radius_m)
        point10 = extend_pnt(pipe_origin_m,vector225,radius_m)
        point11 = extend_pnt(pipe_origin_m,vector315,radius_m)
        point12 = extend_pnt(pipe_end_m,vector45, radius_m)
        point13 = extend_pnt(pipe_end_m,vector135,radius_m)
        point14 = extend_pnt(pipe_end_m,vector225,radius_m)
        point15 = extend_pnt(pipe_end_m,vector315,radius_m)
        point16 = extend_pnt(pipe_mid_m,vector0,  radius_m)
        point17 = extend_pnt(pipe_mid_m,vector90, radius_m)
        point18 = extend_pnt(pipe_mid_m,vector180,radius_m)
        point19 = extend_pnt(pipe_mid_m,vector270,radius_m)
        point20 = extend_pnt(pipe_mid_m,vector315,radius_m)
        point21 = extend_pnt(pipe_mid_m,vector135,radius_m)
        point22 = extend_pnt(pipe_mid_m,vector45, radius_m)
        point23 = extend_pnt(pipe_mid_m,vector225,radius_m)
        pipe_points    = [point00,point01,point02,point03,point04,point05,
                          point06,point07,point08,point09,point10,point11,
                          point12,point13,point14,point15,point16,point17,
                          point18,point19,point20,point21,point22,point23]

        self.vtk_points = pipe_points
        return

    def generate_vtk_points_text(self) -> None:
        """
        Generates the VTK points text.
        """
        text = ""
        i = 0
        for coord in self.vtk_points:
            if (i%3 == 0 ) and (i > 0):
                text += '\r\n'
                new_string = (f"{coord[0]:>.5f} "
                              f"{coord[1]:>.5f} "
                              f"{coord[2]:>.5f} ")
                text += new_string
                i = 1
            else:
                new_string = (f"{coord[0]:>.5f} "
                              f"{coord[1]:>.5f} "
                              f"{coord[2]:>.5f} ")
                text += new_string
                i += 1

        text += '\r\n'

        self.vtk_points_text = text

        return

    def generate_vtk_points_ids_text(self) -> None:
        """
        Generates the VTK points ids text.
        """
        text = ""
        for id_n in self.vtk_point_ids:
            text += f"{id_n:d} "

        text += '\r\n'

        self.vtk_point_ids_text = text

        return

    def is_inside_stl(self, stl_mesh: object) -> bool:
        """
        Checks if the cylinder is inside the STL file.

        Args:
            stl_file (str): The STL file.

        Returns:
            bool: True if the cylinder is inside the STL file, False otherwise.
        """

        check_points = [self.basis_cm, self.end_point_cm]

        points_poly = pv.PolyData(check_points)

        enclosed_points = vtkSelectEnclosedPoints()
        enclosed_points.SetTolerance(0.000001)
        enclosed_points.SetInputData(points_poly)
        enclosed_points.SetSurfaceData(stl_mesh)
        enclosed_points.Update()
        inside_filter = enclosed_points.GetOutput().GetPointData().GetArray("SelectedPoints")
        inside_filter = np.array([inside_filter.GetValue(i) for i in range(inside_filter.GetNumberOfTuples())])
        if not np.array_equal(inside_filter, [0, 0]):
            return True
        else:
            return False

    def __str__(self) -> str:
        return (
            f"CDGSCylinderNode object \n"
            f"PROPERTIES:\n"
            f"mesh id {self.mesh_id}\n"
            f"intensity {self.node_intensity_s:.2e}\n"
            f"basis {self.basis_cm}\n"
            f"end point {self.end_point_cm}\n"
            f"axis 1 {self.axis_1}\n"
            f"radius {self.radius_cm}\n"
            f"height {self.height_cm}\n"
            f"volume {self.node_volume_cm3:.2e}\n"
            f"lateral surface area {self.lateral_surface_area_cm2:.2e}\n"
        )
