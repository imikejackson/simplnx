import simplnx as nx

import itkimageprocessing as cxitk
import orientationanalysis as cxor
import simplnx_test_dirs as nxtest

import numpy as np

# Create a Data Structure
data_structure = nx.DataStructure()

# Filter 1
# Instantiate Filter
nx_filter = nx.ReadStlFileFilter()
# Execute Filter with Parameters
result = nx_filter.execute(
    data_structure=data_structure,
    face_attribute_matrix_name="Face Data",
    face_normals_name="Face Normals",
    scale_output=False,
    stl_file_path=nxtest.get_data_directory() / "STL_Models/ASTMD638_specimen.stl",
    output_triangle_geometry_path=nx.DataPath("Blade"),
    vertex_attribute_matrix_name="Vertex Data"
)
nxtest.check_filter_result(nx_filter, result)

# Filter 2
# Instantiate Filter
nx_filter = nx.ReadStlFileFilter()
# Execute Filter with Parameters
result = nx_filter.execute(
    data_structure=data_structure,
    face_attribute_matrix_name="Face Data",
    face_normals_name="Face Normals",
    #scale_factor: float = ...,
    scale_output=False,
    stl_file_path=nxtest.get_data_directory() / "STL_Models/ASTMD638_specimen.stl",
    output_triangle_geometry_path=nx.DataPath("Blade [Rotation]"),
    vertex_attribute_matrix_name="Vertex Data"
)
nxtest.check_filter_result(nx_filter, result)

# Filter 3
# Instantiate Filter
nx_filter = nx.ApplyTransformationToGeometryFilter()
# Execute Filter with Parameters
result = nx_filter.execute(
    data_structure=data_structure,
    #cell_attribute_matrix_path: DataPath = ...,
    #computed_transformation_matrix: DataPath = ...,
    interpolation_type_index=2,
    #manual_transformation_matrix: List[List[float]] = ...,
    rotation=[1.0, 0.0, 0.0, 180],
    #scale: List[float] = ...,
    input_image_geometry_path=nx.DataPath("Blade [Rotation]"),
    transformation_type_index=3,
    translate_geometry_to_global_origin=False
    #translation: List[float] = ...
)
nxtest.check_filter_result(nx_filter, result)

# Filter 4
# Instantiate Filter
nx_filter = nx.ReadStlFileFilter()
# Execute Filter with Parameters
result = nx_filter.execute(
    data_structure=data_structure,
    face_attribute_matrix_name="Face Data",
    face_normals_name="Face Normals",
    scale_output=False,
    stl_file_path=nxtest.get_data_directory() / "STL_Models/ASTMD638_specimen.stl",
    output_triangle_geometry_path=nx.DataPath("Blade [Translation]"),
    vertex_attribute_matrix_name="Vertex Data"
)
nxtest.check_filter_result(nx_filter, result)

# Filter 5
# Instantiate Filter
nx_filter = nx.ApplyTransformationToGeometryFilter()
# Execute Filter with Parameters
result = nx_filter.execute(
    data_structure=data_structure,
    interpolation_type_index=2,
    input_image_geometry_path=nx.DataPath("Blade [Translation]"),
    transformation_type_index=4,
    translate_geometry_to_global_origin=False,
    translation=[10.0, 10.0, 10.0]
)
nxtest.check_filter_result(nx_filter, result)

# Filter 6
# Instantiate Filter
nx_filter = nx.ReadStlFileFilter()
# Execute Filter with Parameters
result = nx_filter.execute(
    data_structure=data_structure,
    face_attribute_matrix_name="Face Data",
    face_normals_name="Face Normals",
    scale_output=False,
    stl_file_path=nxtest.get_data_directory() / "STL_Models/ASTMD638_specimen.stl",
    output_triangle_geometry_path=nx.DataPath("Blade [Scale]"),
    vertex_attribute_matrix_name="Vertex Data"
)
nxtest.check_filter_result(nx_filter, result)

# Filter 7
# Instantiate Filter
nx_filter = nx.ApplyTransformationToGeometryFilter()
# Execute Filter with Parameters
result = nx_filter.execute(
    data_structure=data_structure,
    interpolation_type_index=2,
    scale=[0.5, 0.5, 1.0],
    input_image_geometry_path=nx.DataPath("Blade [Scale]"),
    transformation_type_index=5,
    translate_geometry_to_global_origin=False
)
nxtest.check_filter_result(nx_filter, result)

# Filter 8
# Instantiate Filter
nx_filter = nx.WriteDREAM3DFilter()
# Set Output File Path
output_file_path = nxtest.get_data_directory() / "Output/ApplyTransformation_Node.dream3d"
# Execute Filter with Parameters
result = nx_filter.execute(
    data_structure=data_structure,
    export_file_path=output_file_path,
    write_xdmf_file=True
)
nxtest.check_filter_result(nx_filter, result)

# *****************************************************************************
# THIS SECTION IS ONLY HERE FOR CLEANING UP THE CI Machines
# If you are using this code, you should COMMENT out the next line
nxtest.cleanup_test_file(output_file_path)
# *****************************************************************************


print("===> Pipeline Complete")
