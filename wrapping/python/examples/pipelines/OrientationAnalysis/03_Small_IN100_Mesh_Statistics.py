import simplnx as nx

import itkimageprocessing as cxitk
import orientationanalysis as cxor
import simplnx_test_dirs as nxtest

import numpy as np

# Create a Data Structure
data_structure = nx.DataStructure()

# Filter 1
# Instantiate Import Data Parameter
import_data = nx.Dream3dImportParameter.ImportData()
import_data.file_path = str(nxtest.get_data_directory() / "Output/SurfaceMesh/SmallIN100_Smoothed.dream3d")
import_data.data_paths = None
# Instantiate Filter
nx_filter = nx.ReadDREAM3DFilter()
# Execute Filter with Parameters
result = nx_filter.execute(data_structure=data_structure, import_data_object=import_data)
nxtest.check_filter_result(nx_filter, result)

# Filter 2
# Instantiate Filter
nx_filter = nx.ComputeTriangleAreasFilter()
# Execute Filter with Parameters
result = nx_filter.execute(
    data_structure=data_structure,
    triangle_areas_array_name="FaceAreas",
    input_triangle_geometry_path=nx.DataPath("TriangleDataContainer")
)
nxtest.check_filter_result(nx_filter, result)

# Filter 3
# Instantiate Filter
nx_filter = nx.TriangleNormalFilter()
# Execute Filter with Parameters
result = nx_filter.execute(
    data_structure=data_structure,
    output_normals_array_name="FaceNormals",
    input_triangle_geometry_path=nx.DataPath("TriangleDataContainer")
)
nxtest.check_filter_result(nx_filter, result)

# Filter 4
# Instantiate Filter
nx_filter = nx.TriangleDihedralAngleFilter()
# Execute Filter with Parameters
result = nx_filter.execute(
    data_structure=data_structure,
    surface_mesh_triangle_dihedral_angles_array_name="FaceDihedralAngles",
    input_triangle_geometry_path=nx.DataPath("TriangleDataContainer")
)
nxtest.check_filter_result(nx_filter, result)

# Filter 5
# Instantiate Filter
nx_filter = cxor.ComputeFaceIPFColoringFilter()
# Execute Filter with Parameters
result = nx_filter.execute(
    data_structure=data_structure,
    crystal_structures_array_path=nx.DataPath("DataContainer/Cell Ensemble Data/CrystalStructures"),
    feature_euler_angles_array_path=nx.DataPath("DataContainer/Cell Feature Data/AvgEulerAngles"),
    feature_phases_array_path=nx.DataPath("DataContainer/Cell Feature Data/Phases"),
    surface_mesh_face_ipf_colors_array_name="FaceIPFColors",
    surface_mesh_face_labels_array_path=nx.DataPath("TriangleDataContainer/Face Data/FaceLabels"),
    surface_mesh_face_normals_array_path=nx.DataPath("TriangleDataContainer/Face Data/FaceNormals")
)
nxtest.check_filter_result(nx_filter, result)

# Filter 6
# Instantiate Filter
nx_filter = cxor.ComputeFeatureFaceMisorientationFilter()
# Execute Filter with Parameters
result = nx_filter.execute(
    data_structure=data_structure,
    avg_quats_array_path=nx.DataPath("DataContainer/Cell Feature Data/AvgQuats"),
    crystal_structures_array_path=nx.DataPath("DataContainer/Cell Ensemble Data/CrystalStructures"),
    feature_phases_array_path=nx.DataPath("DataContainer/Cell Feature Data/Phases"),
    surface_mesh_face_labels_array_path=nx.DataPath("TriangleDataContainer/Face Data/FaceLabels"),
    surface_mesh_face_misorientation_colors_array_name="FaceMisorientationColors"
)
nxtest.check_filter_result(nx_filter, result)

# Filter 7
# Instantiate Filter
nx_filter = nx.WriteDREAM3DFilter()
# Execute Filter with Parameters
output_file_path = nxtest.get_data_directory() / "Output/SurfaceMesh/SmallIN100_MeshStats.dream3d"
result = nx_filter.execute(
    data_structure=data_structure,
    export_file_path=output_file_path,
    write_xdmf_file=True
)
# Check result for errors or warnings
nxtest.check_filter_result(nx_filter, result)

# *****************************************************************************
# THIS SECTION IS ONLY HERE FOR CLEANING UP THE CI Machines
# If you are using this code, you should COMMENT out the next line
nxtest.cleanup_test_file(import_data.file_path)
# *****************************************************************************

print("===> Pipeline Complete")
