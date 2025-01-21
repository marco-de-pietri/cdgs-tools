"""
Module for the CDGS node class.
"""
import math
import pyvista as pv
from vtk import vtkSelectEnclosedPoints
import numpy as np

class CDGSNode:
    """
    class for cdgs node
    """

    class CDGSNode:
        """
        CDGS node generic class.
        """
        def __init__(self) -> None:
            self._mesh_id = 0
            self._energy_bins = []
            self._energy_bin_vals = []
            self._node_intensity = 0
            self._node_comment = ""
            self._cooling_time = 0
            self._cdgs_text = ""

        @property
        def mesh_id(self) -> int:
            """
            Returns the mesh ID.
            """
            return self._mesh_id

        @mesh_id.setter
        def mesh_id(self, value: int) -> None:
            self._mesh_id = value

        @property
        def cooling_time(self) -> float:
            """
            Returns the cooling time.
            """
            return self._cooling_time

        @cooling_time.setter
        def cooling_time(self, value: float) -> None:
            self._cooling_time = value

        @property
        def node_comment(self) -> str:
            """
            Returns the string comment.
            """
            return self._node_comment

        @node_comment.setter
        def node_comment(self, value: str) -> None:
            self._node_comment = value

        @property
        def node_intensity(self) -> float:
            """
            Returns the intensity of the node.
            """
            return self._node_intensity

        @node_intensity.setter
        def node_intensity(self, value: float) -> None:
            self._node_intensity = value

        @property
        def cdgs_text(self) -> str:
            """
            Returns the CDGS text.
            """
            return self._cdgs_text

        @cdgs_text.setter
        def cdgs_text(self, value: str) -> None:
            """
            set the cdgs text
            """
            self._cdgs_text = value



        def __str__(self) -> str:
            return (
                f"CDGSNode object with mesh id {self.mesh_id} "
                f"and node intensity {self.node_intensity:.2e}"
            )
class CDGSCylinderNode(CDGSNode):
    """
    class for cdgs cylinder node
    """

    def __init__(self) -> None:
        super().__init__()
        self._basis = []
        self._end_point = []
        self._axis = []
        self._radius = 0
        self._height = 0
        self._node_volume = 0
        self._lateral_surface_area = 0

    @property
    def basis(self) -> list:
        """
        Returns the basis list.

        :return: A list representing the basis.
        :rtype: list
        """
        return self._basis

    @basis.setter
    def basis(self, value: list) -> None:
        self._basis = value

    @property
    def end_point(self) -> list:
        """
        Returns the end point list.

        :return: A list representing the end point.
        :rtype: list
        """
        return self._end_point

    def calculate_end_point(self) -> None:
        """
        Calculates the end point of the cylinder.
        """
        self._end_point = [
            self.basis[0] + self.axis[0] * self.height,
            self.basis[1] + self.axis[1] * self.height,
            self.basis[2] + self.axis[2] * self.height,
        ]

    def generate_cdgs_text(self) -> None:
        """
        Generates the CDGS text.
        """
        text = ""
        text += f"mesh_id {self.mesh_id}\n"
        text += f"{self.node_comment}\n"
        text += f"Cooling_time {self.cooling_time:.6e}\n"
        text += f"total_source {self.node_intensity:.9e}\n"
        text += f"energy_type bins\n"
        text += f"energy_boundaries {len(self.energy_bins):d}\n"

        self.cdgs_text = text

        return

    @property
    def axis(self) -> list:
        """
        Returns the axis attribute.

        Returns:
            list: The axis attribute.
        """
        return self._axis

    @axis.setter
    def axis(self, value: list) -> None:
        self._axis = value

    @property
    def radius(self) -> float:
        """
        Returns the radius of the node.

        Returns:
            float: The radius of the node.
        """
        return self._radius

    @radius.setter
    def radius(self, value: float) -> None:
        self._radius = value

    @property
    def height(self) -> float:
        """
        Returns the length of the object.

        Returns:
            float: The length of the object.
        """
        return self._height

    @height.setter
    def height(self, value: float) -> None:
        self._height = value

    @property
    def node_volume(self) -> float:
        """
        Returns the volume of the cylinder.

        Returns:
            float: The volume of the cylinder.
        """
        return self._node_volume

    def calculate_volume(self) -> None:
        """
        Calculates the volume of the cylinder.
        """
        self._node_volume = math.pi * self.height * (self.radius ** 2)

    @property
    def lateral_surface_area(self) -> float:
        """
        Returns the lateral surface area of the cylinder.

        Returns:
            float: The lateral surface area of the cylinder.
        """
        return self._lateral_surface_area

    def calculate_lateral_surface_area(self) -> None:
        """
        Calculates the lateral surface area of the cylinder.
        """
        self._lateral_surface_area = 2 * math.pi * self.radius * self.height

    def is_inside_stl(self, stl_mesh: object) -> bool:
        """
        Checks if the cylinder is inside the STL file.

        Args:
            stl_file (str): The STL file.

        Returns:
            bool: True if the cylinder is inside the STL file, False otherwise.
        """

        check_points = [self.basis, self.end_point]

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
            f"intensity {self.node_intensity:.2e}\n"
            f"basis {self.basis}\n"
            f"end point {self.end_point}\n"
            f"axis {self.axis}\n"
            f"radius {self.radius}\n"
            f"height {self.height}\n"
            f"volume {self.node_volume:.2e}\n"
            f"lateral surface area {self.lateral_surface_area:.2e}\n"
        )
