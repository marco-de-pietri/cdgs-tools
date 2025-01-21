"""
module for cdgs class
"""

import copy
import re
import math
import sys
import pyvista as pv

from .cdgs_node_class import CDGSCylinderNode

def get_node_comment(text):
    """
    function that parse the cdgs text to extract the node comment
    """

    commentPat = re.compile (r"^mesh_id.{1,}?(?=^Cooling_time)",re.IGNORECASE | re.MULTILINE | re.DOTALL)
    commentLine = commentPat.findall(text)
    return commentLine[0].splitlines()[1].strip(" \r\n")

def get_cooling_time(text):
    """
    function that parse the cdgs text to extract the node cooling time
    """


    cooling_time_pat = re.compile(r"[\r\n]cooling_time\s+.*", re.IGNORECASE )


    ct_matches = cooling_time_pat.findall(text)


    return float(ct_matches[0].split()[1])


def get_cylinder_geometry(text):
    """
    this function parse the text of the cdgs file to extract the node geometry
    at the moment it assumes we are working with a cylinder geometry
    """

    geometry_pat = re.compile(
        r"^mesh_type.{1,}?(?=^source_data)", re.IGNORECASE | re.MULTILINE | re.DOTALL
    )
    check_geo_1_pat = re.compile(r"^mesh_type\s+.*", re.IGNORECASE | re.MULTILINE)
    check_geo_2_pat = re.compile(r"^mesh_boundaries\s+.*", re.IGNORECASE | re.MULTILINE)

    geom_text = geometry_pat.findall(text)
    check_geo_1_line = check_geo_1_pat.findall(geom_text[0])
    check_geo_2_line = check_geo_2_pat.findall(geom_text[0])

    if check_geo_1_line[0].split()[1] != "cyl":
        print("error mesh is not cyl!")
        sys.exit()

    if (
        (check_geo_2_line[0].split()[1] != "2")
        and (check_geo_2_line[0].split()[2] != "2")
        and (check_geo_2_line[0].split()[3] != "2")
    ):
        print("error with mesh boundaries")
        sys.exit()

    geometry_lines = geom_text[0].splitlines()

    basis = [float(value) for value in geometry_lines[2].split()]
    axis_1 = [float(value) for value in geometry_lines[3].split()]
    axis_2 = [float(value) for value in geometry_lines[4].split()]

    # norm = pow(sum([pow(val,2) for val in axis]),0.5)

    radius = float(geometry_lines[5].split()[1])
    length = float(geometry_lines[7].split()[1])

    return basis, axis_1, axis_2, radius, length


def get_node_intensity(text):
    """
    function that parse the cdgs text to extract the intensity of the cdgs node
    """

    pat_intensity = re.compile(r"^total_source\s+.*", re.IGNORECASE | re.MULTILINE)

    intensity_line = pat_intensity.findall(text)

    return float(intensity_line[0].split()[1])


def get_mesh_id(text):
    """
    function that parse the cdgs text to extract the mesh id
    """

    pat_mesh_id = re.compile(r"^mesh_id\s+.*", re.IGNORECASE )

    mesh_line = pat_mesh_id.findall(text)

    return int(mesh_line[0].split()[1])


def get_energy_bin_vals(text):
    """
    function that parse the cdgs text to extract the energy bin values
    """
    pat_ebin_vals = re.compile(
        r"^source_data.{1,}?(?=\Z)", re.IGNORECASE | re.MULTILINE | re.DOTALL
    )

    ebin_val_lines = pat_ebin_vals.findall(text)[0].splitlines()[3:]

    ebin_val_list = []

    for line in ebin_val_lines:
        values = line.split()
        for val in values:
            ebin_val_list.append(float(val))

    mid_index = len(ebin_val_list) // 2
    ebin_intensity_list = ebin_val_list[:mid_index]
    ebin_error_list = ebin_val_list[mid_index:]

    return ebin_intensity_list, ebin_error_list


def get_energy_bins(text):
    """
    function that parse the cdgs text to extract the energy boundaries
    """

    pat_n_ebins = re.compile(r"^energy_boundaries\s+.*", re.IGNORECASE | re.MULTILINE)

    pat_ebins = re.compile(
        r"^energy_boundaries.{1,}?(?=^mesh_type)",
        re.IGNORECASE | re.MULTILINE | re.DOTALL,
    )

    n_ebins_line = pat_n_ebins.findall(text)
    ebins_lines = pat_ebins.findall(text)[0].splitlines()[1:]
    eb_list = []
    for line in ebins_lines:
        values = line.split()
        for val in values:
            eb_list.append(float(val))

    n_ebins = int(n_ebins_line[0].split()[1])

    if len(eb_list) != n_ebins:
        print("mismatch number of energy boundaries!")
        sys.exit()

    return eb_list


def get_file_tot_intesity(path):
    """
    parse the text file to get the total intensity listed in the file header
    """

    rf = open(path, "r", encoding="utf8", errors="ignore")
    for _, line in enumerate(rf):
        if line.split()[0] == "global_source":
            tot_intensity = float(line.split()[1])
            break

    return tot_intensity


def get_file_mesh_n(path):
    """
    parse the text file to get the number of cdgs nodes
    """
    rf = open(path, "r", encoding="utf8", errors="ignore")
    for _, line in enumerate(rf):
        if line.split()[0] == "num_meshes":
            tot_mesh = int(line.split()[1])
            break
    return tot_mesh


def read_cdgs_file(cdgs_file: str) -> object:
    """
    read the cdgs file and return a CDGS object
    """
    block_pat = re.compile(
        "(?=^mesh_id).{1,}?(?=^mesh_id|^end_source_data)",
        re.IGNORECASE | re.MULTILINE | re.DOTALL,
    )

    cdgs_nodes_number = get_file_mesh_n(cdgs_file)
    mesh_tot_intensity = get_file_tot_intesity(cdgs_file)


    try:
        imp_file = open(cdgs_file, "r", encoding="utf8", errors="ignore")
    except IOError:
        print("couldn't open DGS file")
        sys.exit()
    with imp_file:

        cdgs_object = CDGS()

        text = imp_file.read()
        node_blocks = block_pat.findall(text)

        for node_text in node_blocks:

            cyl_node = CDGSCylinderNode()
            mesh_id = get_mesh_id(node_text)
            cyl_node.mesh_id = mesh_id

            node_intensity = get_node_intensity(node_text)
            cyl_node.node_intensity_s = node_intensity

            node_comment = get_node_comment(node_text)
            cyl_node.node_comment = node_comment

            node_cooling_time = get_cooling_time(node_text)
            cyl_node.cooling_time = node_cooling_time

            energy_bins = get_energy_bins(node_text)
            cyl_node.energy_bins = energy_bins

            energy_bin_vals, energy_bin_errs = get_energy_bin_vals(node_text)
            cyl_node.energy_bin_vals = energy_bin_vals
            cyl_node.energy_bin_errs = energy_bin_errs


            basis, axis_1, axis_2, radius, height = get_cylinder_geometry(node_text)
            cyl_node.basis_cm = basis
            cyl_node.axis_1 = axis_1
            cyl_node.axis_2 = axis_2
            cyl_node.radius_cm = radius
            cyl_node.height_cm = height

            cyl_node.calculate_volume()
            cyl_node.calculate_lateral_surface_area()
            cyl_node.calculate_end_point()
            cyl_node.generate_cdgs_text()
            cyl_node.generate_vtk_points()
            cyl_node.generate_vtk_points_text()
            point_ids = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
                          12,13,14,15,16,17,18,19,20,21,22,23]
            cyl_node.vtk_point_ids = [(cdgs_object._tot_vtk_points + i) for i in point_ids]
            cyl_node.generate_vtk_points_ids_text()
            cyl_node.vtk_cell_type = 33
            cyl_node.node_intensity_m3_s = cyl_node.node_intensity_s / cyl_node.node_volume_m3

            cdgs_object.add_node(cyl_node)


        if cdgs_object.tot_meshes != cdgs_nodes_number:
            print("mismatch number of nodes!")
            print ("nodes parsed: ", cdgs_object.tot_meshes)
            print ("nodes in file: ", cdgs_nodes_number)
            sys.exit()
        if not math.isclose(cdgs_object.tot_intensity, mesh_tot_intensity, rel_tol=1e-3):
            print("mismatch total intensity!")
            print(f"intensity parsed: {cdgs_object.tot_intensity:.8e}")
            print (f"intensity in file header: , {mesh_tot_intensity:.8e}")
            sys.exit()

    return cdgs_object


class CDGS:
    """
    class for cdgs
    """

    def __init__(self) -> None:
        self._tot_meshes = 0
        self._tot_intensity = 0
        self._tot_volume = 0
        self._nodes = []
        self._tot_vtk_points = 0

    def __str__(self) -> str:
        return (
            f"CDGS object\n"
            f"PROPERTY: \n"
            f"total meshes {self.tot_meshes} \n"
            f"total intensity [Bq] {self.tot_intensity:.2e}"
            f"total volume [cm3] {self.tot_volume:.2e}"
        )

    @property
    def nodes(self) -> list:
        """
        Returns the list of nodes.
        """
        return self._nodes

    @property
    def tot_meshes(self) -> int:
        """
        Returns the total number of meshes.
        """
        return self._tot_meshes

    @property
    def tot_intensity(self) -> float:
        """
        Returns the total intensity.
        """
        return self._tot_intensity

    @property
    def tot_vtk_points(self) -> float:
        """
        number of vtk points
        """
        return self._tot_vtk_points

    @property
    def tot_volume(self) -> float:
        """
        Returns the total volume
        """
        return self._tot_volume


    def add_node(self, node) -> None:
        """
        Add a node to the list of nodes.
        """
        self._nodes.append(copy.deepcopy(node))
        self._tot_meshes += 1
        self._tot_intensity += node.node_intensity_s
        self._tot_volume += node.node_volume_cm3
        self._tot_vtk_points += len(node.vtk_points)

    def renumber_nodes(self) -> None:
        """
        renumber the nodes in the cdgs object
        """
        for i, node in enumerate(self.nodes):
            node.mesh_id = i

    def filter_by_stl (self, stl_file) -> object:
        """
        filter the cdgs object leaving only the nodes that are inside the stl file
        """

        mesh = pv.read(stl_file)
        cdgs_subset = CDGS()

        for node in self.nodes:

            if node.is_inside_stl( mesh):
                cdgs_subset.add_node(node)

        cdgs_subset.renumber_nodes()

        return cdgs_subset

    def write_cdgs (self, output_file) -> None:
        """
        write the cdgs object to a file
        """

        with open(output_file, "w") as wf:
            wf.write(f"num_meshes {self.tot_meshes}\n")
            wf.write(f"global_source {self.tot_intensity:.7e}\n")

            for node in self.nodes:
                wf.write(node.cdgs_text)

    def write_vtk(self,file_name):
        """
        this function writes the vtk file of the pipe nodes for visualization
        """

        n_vtk_points = sum([len(node.vtk_points) for node in self.nodes])

        n_nodes = len(self.nodes)

        try:
            out = open(file_name,'w',encoding="utf8", errors='ignore')
        except IOError:
            raise IOError("couldn't write vtk file")
        with out:

            header = """# vtk DataFile Version 2.0
DGS printing
ASCII
DATASET UNSTRUCTURED_GRID\r\n"""
            out.write(header)
            out.write(f"POINTS  {n_vtk_points:d}  floats\r\n")

            for node in self.nodes:
                out.write(node.vtk_points_text)

            out.write(f"CELLS  {n_nodes:d}   {(n_nodes + n_vtk_points):d}\r\n")

            for node in self.nodes:
                n_point_ids = len(node.vtk_points)
                intro_string = f"{n_point_ids:d} "
                write_string =  intro_string + node.vtk_point_ids_text
                out.write(write_string)

            out.write(f"CELL_TYPES  {n_nodes:d}  \r\n")

            for node in self.nodes:
                out.write(f"{node.vtk_cell_type:d}\r\n")

            out.write(f"CELL_DATA  {n_nodes:d} \r\n")

            out.write("SCALARS node_id int 1 \r\n")
            out.write("LOOKUP_TABLE default\r\n")

            for node in self.nodes:
                out.write(f"{node.mesh_id:d}\r\n")

            out.write("SCALARS node_emission_intensity_particle_s double 1 \r\n")
            out.write("LOOKUP_TABLE default\r\n")
            for node in self.nodes:
                out.write(f"{node.node_intensity_s:.5e}\r\n")

            out.write("SCALARS node_emission_rate_particle_m3_s double 1 \r\n")
            out.write("LOOKUP_TABLE default\r\n")
            for node in self.nodes:
                out.write(f"{node.node_intensity_m3_s:.5e}\r\n")

            out.write("SCALARS volume_m3 double 1 \r\n")
            out.write("LOOKUP_TABLE default\r\n")
            for node in self.nodes:
                out.write(f"{node.node_volume_m3:.5e}\r\n")

        return





