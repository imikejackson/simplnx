# https://www.hmilch.net/downloads/cli_format.html

import simplnx as nx
import numpy as np
from enum import Enum
import sys
import time
import copy
import re
from typing import List, Tuple
from pathlib import Path
from .common.Result import Result, make_error_result

class GeometricFeature(object):
  def __init__(self, feature_id, layer_id, z_height, data: dict, poly_id, n, start_xvals, start_yvals, end_xvals, end_yvals, dir = None) -> Result:
        self.feature_id = feature_id
        self.layer_id = layer_id
        self.z_height = z_height
        self.poly_id = poly_id
        self.dir = dir
        self.n = n
        self.start_xvals = start_xvals
        self.start_yvals = start_yvals
        self.end_xvals   = end_xvals
        self.end_yvals   = end_yvals
        self.data = data

class CliReaderFilter:
  # Parameter Keys
  CLI_FILE_PATH_KEY = 'cli_file_path'
  USE_X_DIMENSION_RANGE_KEY = 'use_x_dimension_range'
  USE_Y_DIMENSION_RANGE_KEY = 'use_y_dimension_range'
  USE_Z_DIMENSION_RANGE_KEY = 'use_z_dimension_range'
  MIN_MAX_X_COORDS_KEY = 'min_max_x_coords'
  MIN_MAX_Y_COORDS_KEY = 'min_max_y_coords'
  MIN_MAX_Z_COORDS_KEY = 'min_max_z_coords'
  OUT_OF_BOUNDS_BEHAVIOR_KEY = 'out_of_bounds_behavior'
  READ_EXTRA_METADATA = 'read_extra_metadata'
  CREATE_EDGE_TYPE_ARRAY = 'create_geometric_labels_array'
  OUTPUT_EDGE_GEOM_PATH_KEY = 'output_edge_geom_path'
  OUTPUT_VERTEX_ATTRMAT_NAME_KEY = 'output_vertex_attrmat_name'
  OUTPUT_EDGE_ATTRMAT_NAME_KEY = 'output_edge_attrmat_name'
  OUTPUT_FEATURE_ATTRMAT_NAME_KEY = 'output_feature_attrmat_name'
  SHARED_VERTICES_ARRAY_NAME_KEY = 'shared_vertices_array_name'
  SHARED_EDGES_ARRAY_NAME_KEY = 'shared_edges_array_name'

  # Constants
  LAYER_ARRAY_NAME = 'Layer'
  LABEL_ARRAY_NAME = 'Label'
  EDGE_TYPE_ARRAY_NAME = 'Edge Types'

  class OutOfBoundsBehavior(Enum):
    InterpolateOutsideVertex = 0
    IgnoreEdge = 1
    FilterError = 2

  # Out-Of-Bounds Choices
  OUT_OF_BOUNDS_BEHAVIOR_CHOICES = [
    "Interpolate Outside Vertex",
    "Ignore Edge",
    "Filter Error"
  ]

  def uuid(self) -> nx.Uuid:
    return nx.Uuid('ef2e9ad6-862c-4142-982f-704618a8855c')

  def human_name(self) -> str:
    return 'Read CLI File (Python)'

  def class_name(self) -> str:
    return 'CliReaderFilter'

  def name(self) -> str:
    return 'CliReaderFilter'

  def default_tags(self) -> List[str]:
    return ['python']

  def clone(self):
    return CliReaderFilter()

  def parameters(self) -> nx.Parameters:
    params = nx.Parameters()

    params.insert(nx.Parameters.Separator("Parameters"))
    params.insert(nx.FileSystemPathParameter(CliReaderFilter.CLI_FILE_PATH_KEY, 'Input CLI File', 'The path to the input CLI file that will be read.', '', {'.cli'}, nx.FileSystemPathParameter.PathType.InputFile))
    params.insert_linkable_parameter(nx.BoolParameter(CliReaderFilter.USE_X_DIMENSION_RANGE_KEY, 'Use X Dimension Range', 'Determines whether or not to use X bounds to import a range of data in the X dimension.', False))
    params.insert(nx.VectorFloat64Parameter(CliReaderFilter.MIN_MAX_X_COORDS_KEY, 'X Min/Max', 'The minimum and maximum X coordinate for the X bounds.', [0.0, 100.0], ['X Min', 'X Max']))
    params.insert_linkable_parameter(nx.BoolParameter(CliReaderFilter.USE_Y_DIMENSION_RANGE_KEY, 'Use Y Dimension Range', 'Determines whether or not to use Y bounds to import a range of data in the Y dimension.', False))
    params.insert(nx.VectorFloat64Parameter(CliReaderFilter.MIN_MAX_Y_COORDS_KEY, 'Y Min/Max', 'The minimum and maximum Y coordinate for the Y bounds.', [0.0, 100.0], ['Y Min', 'Y Max']))
    params.insert_linkable_parameter(nx.BoolParameter(CliReaderFilter.USE_Z_DIMENSION_RANGE_KEY, 'Use Z Dimension Range', 'Determines whether or not to use Z bounds to import a range of data in the Z dimension.', False))
    params.insert(nx.VectorFloat64Parameter(CliReaderFilter.MIN_MAX_Z_COORDS_KEY, 'Z Min/Max', 'The minimum and maximum Z coordinate for the Z bounds.', [0.0, 100.0], ['Z Min', 'Z Max']))
    params.insert(nx.ChoicesParameter(CliReaderFilter.OUT_OF_BOUNDS_BEHAVIOR_KEY, 'Out-Of-Bounds Behavior', 'The behavior to implement if an edge intersects a bound (one vertex is inside, one vertex is outside).\n\n"Interpolate Outside Vertex" will move the outside vertex of a boundary-intersecting edge from its current position to the boundary edge.\n"Ignore Edge" will ignore any edge that intersects a bound.\n"Filter Error" will make this filter throw an error when it encounters an edge that intersects a bound.', 0, CliReaderFilter.OUT_OF_BOUNDS_BEHAVIOR_CHOICES))
    params.insert(nx.BoolParameter(CliReaderFilter.READ_EXTRA_METADATA, 'Read Extra Metadata', 'Determines whether or not to read the extra metadata that may be included (per layer) in the file.', True))
    params.insert(nx.BoolParameter(CliReaderFilter.CREATE_EDGE_TYPE_ARRAY, 'Create Edge Type Array', 'Determines whether or not to create an unsigned 8-bit integer array that labels hatches edges with 0 and polyline edges with 1.', True))
    params.insert(nx.Parameters.Separator("Created Data Objects"))
    params.insert(nx.DataGroupCreationParameter(CliReaderFilter.OUTPUT_EDGE_GEOM_PATH_KEY, 'Output Edge Geometry', 'The path to the newly created edge geometry.', nx.DataPath("[Edge Geometry]")))
    params.insert(nx.DataObjectNameParameter(CliReaderFilter.OUTPUT_VERTEX_ATTRMAT_NAME_KEY, 'Output Vertex Attribute Matrix Name', 'The name of the newly created vertex attribute matrix.', 'Vertex Data'))
    params.insert(nx.DataObjectNameParameter(CliReaderFilter.OUTPUT_EDGE_ATTRMAT_NAME_KEY, 'Output Edge Attribute Matrix Name', 'The name of the newly created edge attribute matrix.', 'Edge Data'))
    params.insert(nx.DataObjectNameParameter(CliReaderFilter.OUTPUT_FEATURE_ATTRMAT_NAME_KEY, 'Output Feature Attribute Matrix Name', 'The name of the newly created feature attribute matrix.', 'Feature Data'))
    params.insert(nx.DataObjectNameParameter(CliReaderFilter.SHARED_VERTICES_ARRAY_NAME_KEY, 'Shared Vertices Array Name', 'The name of the newly created shared vertices array.', 'Shared Vertices'))
    params.insert(nx.DataObjectNameParameter(CliReaderFilter.SHARED_EDGES_ARRAY_NAME_KEY, 'Shared Edges Array Name', 'The name of the newly created shared edges array.', 'Shared Edges'))

    params.link_parameters(CliReaderFilter.USE_X_DIMENSION_RANGE_KEY, CliReaderFilter.MIN_MAX_X_COORDS_KEY, True)
    params.link_parameters(CliReaderFilter.USE_Y_DIMENSION_RANGE_KEY, CliReaderFilter.MIN_MAX_Y_COORDS_KEY, True)
    params.link_parameters(CliReaderFilter.USE_Z_DIMENSION_RANGE_KEY, CliReaderFilter.MIN_MAX_Z_COORDS_KEY, True)
    params.link_parameters(CliReaderFilter.USE_X_DIMENSION_RANGE_KEY, CliReaderFilter.OUT_OF_BOUNDS_BEHAVIOR_KEY, True)
    params.link_parameters(CliReaderFilter.USE_Y_DIMENSION_RANGE_KEY, CliReaderFilter.OUT_OF_BOUNDS_BEHAVIOR_KEY, True)
    params.link_parameters(CliReaderFilter.USE_Z_DIMENSION_RANGE_KEY, CliReaderFilter.OUT_OF_BOUNDS_BEHAVIOR_KEY, True)

    return params

  def parameters_version(self) -> int:
    return 1

  def preflight_impl(self, data_structure: nx.DataStructure, args: dict, message_handler: nx.IFilter.MessageHandler, should_cancel: nx.AtomicBoolProxy) -> nx.IFilter.PreflightResult:
    cli_file_path: str = args[CliReaderFilter.CLI_FILE_PATH_KEY]
    output_edge_geom_path: nx.DataPath = args[CliReaderFilter.OUTPUT_EDGE_GEOM_PATH_KEY]
    output_vertex_attrmat_name: str = args[CliReaderFilter.OUTPUT_VERTEX_ATTRMAT_NAME_KEY]
    output_edge_attrmat_name: str = args[CliReaderFilter.OUTPUT_EDGE_ATTRMAT_NAME_KEY]
    output_feature_attrmat_name: str = args[CliReaderFilter.OUTPUT_FEATURE_ATTRMAT_NAME_KEY]
    shared_vertices_array_name: str = args[CliReaderFilter.SHARED_VERTICES_ARRAY_NAME_KEY]
    shared_edges_array_name: str = args[CliReaderFilter.SHARED_EDGES_ARRAY_NAME_KEY]
    use_x_dimension_range: bool = args[CliReaderFilter.USE_X_DIMENSION_RANGE_KEY]
    use_y_dimension_range: bool = args[CliReaderFilter.USE_Y_DIMENSION_RANGE_KEY]
    use_z_dimension_range: bool = args[CliReaderFilter.USE_Z_DIMENSION_RANGE_KEY]
    min_max_x_coords: list = args[CliReaderFilter.MIN_MAX_X_COORDS_KEY]
    min_max_y_coords: list = args[CliReaderFilter.MIN_MAX_Y_COORDS_KEY]
    min_max_z_coords: list = args[CliReaderFilter.MIN_MAX_Z_COORDS_KEY]
    read_extra_metadata: bool = args[CliReaderFilter.READ_EXTRA_METADATA]
    create_edge_types_array: bool = args[CliReaderFilter.CREATE_EDGE_TYPE_ARRAY]
    
    if use_x_dimension_range and min_max_x_coords[0] > min_max_x_coords[1]:
      return nx.IFilter.PreflightResult(nx.OutputActions(), [nx.Error(-9100, f"Invalid X Dimension Range: The minimum X coordinate ({min_max_x_coords[0]}) is larger than the maximum X coordinate ({min_max_x_coords[1]}).")])
    if use_y_dimension_range and min_max_y_coords[0] > min_max_y_coords[1]:
        return nx.IFilter.PreflightResult(nx.OutputActions(), [nx.Error(-9101, f"Invalid Y Dimension Range: The minimum Y coordinate ({min_max_y_coords[0]}) is larger than the maximum Y coordinate ({min_max_y_coords[1]}).")])
    if use_z_dimension_range and min_max_z_coords[0] > min_max_z_coords[1]:
        return nx.IFilter.PreflightResult(nx.OutputActions(), [nx.Error(-9102, f"Invalid Z Dimension Range: The minimum Z coordinate ({min_max_z_coords[0]}) is larger than the maximum Z coordinate ({min_max_z_coords[1]}).")])

    # Here we create the Edge Geometry (and the 2 internal Attribute Matrix to hold vertex and edge data arrays.)
    # Because this is a "reader" type of filter we do not know (at least in this reader implementation)
    # the number of vertices or edges at preflight time. During execute we will need to ensure that
    # everything is sized correctly.
    output_actions = nx.OutputActions()
    output_actions.append_action(nx.CreateEdgeGeometryAction(geometry_path=output_edge_geom_path, num_edges=1, num_vertices=1, vertex_attribute_matrix_name=output_vertex_attrmat_name, edge_attribute_matrix_name=output_edge_attrmat_name, shared_vertices_name=shared_vertices_array_name, shared_edges_name=shared_edges_array_name))

    array_names, num_of_labels = self._parse_geometry_array_names(Path(cli_file_path))
    # Because extra geometric data is not included in the specification, we are setting 'Layer' array to int32 and all other arrays to float32
    edge_attr_mat_path = output_edge_geom_path.create_child_path(output_edge_attrmat_name)

    preflight_updated_values = []
    if read_extra_metadata:
      preflight_value = nx.IFilter.PreflightValue()
      preflight_value.name = "Extra Metadata Arrays"
      array_names_str = '\n'.join(array_names)
      preflight_value.value = array_names_str
      preflight_updated_values.append(preflight_value)

      for array_name in array_names:
        dtype = nx.DataType.float32
        if array_name == self.LAYER_ARRAY_NAME:
          dtype = nx.DataType.int32
        array_path = edge_attr_mat_path.create_child_path(array_name)
        output_actions.append_action(nx.CreateArrayAction(dtype, [1], [1], array_path))
    
      if num_of_labels > 0:
        array_path = edge_attr_mat_path.create_child_path(self.LABEL_ARRAY_NAME)
        output_actions.append_action(nx.CreateArrayAction(nx.DataType.int32, [1], [1], array_path))

        feature_attr_mat_path = output_edge_geom_path.create_child_path(output_feature_attrmat_name)
        output_actions.append_action(nx.CreateAttributeMatrixAction(feature_attr_mat_path, [num_of_labels]))

        label_array_path = feature_attr_mat_path.create_child_path(self.LABEL_ARRAY_NAME)
        output_actions.append_action(nx.CreateStringArrayAction([num_of_labels], label_array_path))
    
    if create_edge_types_array:
      array_path = edge_attr_mat_path.create_child_path(self.EDGE_TYPE_ARRAY_NAME)
      output_actions.append_action(nx.CreateArrayAction(nx.DataType.uint8, [1], [1], array_path))

    return nx.IFilter.PreflightResult(output_actions, preflight_values=preflight_updated_values)

  def _process_layers(self, layer_features, hatch_labels, message_handler: nx.IFilter.MessageHandler) -> Result:
    last_message_time = 0
    start_vertices = []
    end_vertices = []
    data_arrays = {}
    if len(hatch_labels) > 0:
      data_arrays[self.LABEL_ARRAY_NAME] = []
    data_arrays[self.EDGE_TYPE_ARRAY_NAME] = []
    num_of_edges = 0
    for layer_idx in range(len(layer_features)):
      layer = layer_features[layer_idx]
      if not layer:
        current_time = time.time()
        if message_handler is not None and (current_time - last_message_time) >= 1:
          message_handler(nx.IFilter.Message(nx.IFilter.Message.Type.Info, f'Skipping layer {layer_idx + 1}/{len(layer_features)}...'))
          last_message_time = current_time
          continue

      for feature in layer:
        if feature.n == 0:
          continue

        num_of_edges += feature.n
        if len(feature.start_xvals) != len(feature.start_yvals):
          return Result([nx.Error(-2011, f"The start x value and start y value lists for hatch (id {feature.hatch_id}) are not the same size!")])
        if len(feature.end_xvals) != len(feature.end_yvals):
          return Result([nx.Error(-2012, f"The end x value and end y value lists for hatch (id {feature.hatch_id}) are not the same size!")])
        for start_x, start_y in zip(feature.start_xvals, feature.start_yvals):
          start_vertices.append([start_x, start_y, feature.z_height])
        for end_x, end_y in zip(feature.end_xvals, feature.end_yvals):
          end_vertices.append([end_x, end_y, feature.z_height])
        if len(hatch_labels) > 0:
          data_arrays[self.LABEL_ARRAY_NAME].extend([''] * feature.n)

        for array_name, value in feature.data.items():
          if not array_name in data_arrays:
            data_arrays[array_name] = [value] * feature.n
          else:
            data_arrays[array_name].extend([value] * feature.n)
        
        data_arrays[self.EDGE_TYPE_ARRAY_NAME].extend([feature.feature_id] * feature.n)
      
      current_time = time.time()
      if message_handler is not None and (current_time - last_message_time) >= 1:
        message_handler(nx.IFilter.Message(nx.IFilter.Message.Type.Info, f'Imported layer {layer_idx + 1}/{len(layer_features)}'))
        last_message_time = current_time

    message_handler(nx.IFilter.Message(nx.IFilter.Message.Type.Info, f'Imported layer {len(layer_features)}/{len(layer_features)}'))
    return Result(value=(start_vertices, end_vertices, data_arrays, num_of_edges))
  
  def _update_edge_geometry(self, data_structure: nx.DataStructure, output_edge_geom_path: nx.DataPath, output_edge_attrmat_name: str,
                            output_feature_attrmat_name:str, start_vertices: list, end_vertices: list, data_arrays: list, num_of_edges: int,
                            hatch_labels: list, read_extra_metadata: bool, create_edge_types_array: bool, message_handler: nx.IFilter.MessageHandler):
    # Update the edge geometry using the start and end vertices, data arrays, and number of edges
    edge_geom: nx.EdgeGeom = data_structure[output_edge_geom_path]

    # Tell the Edge Geometry to resize the shared vertex list so that we can 
    # copy in the vertices.
    message_handler(nx.IFilter.Message(nx.IFilter.Message.Type.Info, f'Saving Vertex List...'))
    vertex_list = [item for pair in zip(start_vertices, end_vertices) for item in pair]
    edge_geom.resize_vertices(len(vertex_list))

    if (len(vertex_list) > 0):
      vertices_array = edge_geom.vertices
      vertices_view = vertices_array.store.npview()
      vertices_view[:] = vertex_list

    # Tell the Edge Geometry to resize the shared edge list so that we can
    # copy in the edge list and also copy in all the edge arrays
    message_handler(nx.IFilter.Message(nx.IFilter.Message.Type.Info, f'Saving Edges...'))
    edge_geom.resize_edges(num_of_edges)

    if (num_of_edges > 0 and len(vertex_list) > 0):
      edges_array = edge_geom.edges
      edges_view = edges_array.store.npview()
      edges_view[:] = [[i, i+1] for i in range(0, len(vertex_list), 2)]

    # Get the nx.DataPath to the Edge Attribute Matrix
    edge_attr_mat_path = output_edge_geom_path.create_child_path(output_edge_attrmat_name)

    if read_extra_metadata:
      # Copy the all the edge data into the edge attribute matrix
      for array_name, values in data_arrays.items():
        if array_name == self.LABEL_ARRAY_NAME:
          continue

        message_handler(nx.IFilter.Message(nx.IFilter.Message.Type.Info, f"Saving Cell Array '{array_name}'..."))
        array_path = edge_attr_mat_path.create_child_path(array_name)
        array: nx.IDataArray = data_structure[array_path]
        values_arr = np.array(values)
        values_arr = values_arr.reshape([len(values)] + array.cdims)

        if len(values) > 0:
          array_view = array.store.npview()
          values_arr = values_arr.astype(float).astype(array_view.dtype)
          array_view[:] = values_arr

      # Save the feature level data
      if len(hatch_labels) > 0:
        feature_attr_mat_path = output_edge_geom_path.create_child_path(output_feature_attrmat_name)
        label_feature_array_path = feature_attr_mat_path.create_child_path(self.LABEL_ARRAY_NAME)
        label_feature_array: nx.StringArray = data_structure[label_feature_array_path]
        message_handler(nx.IFilter.Message(nx.IFilter.Message.Type.Info, f"Saving Feature Array '{self.LABEL_ARRAY_NAME}'..."))
        label_feature_array.initialize_with_list(list(hatch_labels.values()))
    
    if create_edge_types_array:
      message_handler(nx.IFilter.Message(nx.IFilter.Message.Type.Info, f"Saving Edge Type Array..."))
      array_path = edge_attr_mat_path.create_child_path(self.EDGE_TYPE_ARRAY_NAME)
      array: nx.IDataArray = data_structure[array_path]
      values = data_arrays[self.EDGE_TYPE_ARRAY_NAME]
      values_arr = np.array(values)
      values_arr = values_arr.reshape([len(values)] + array.cdims)
      if len(values) > 0:
        array_view = array.store.npview()
        values_arr = values_arr.astype(float).astype(array_view.dtype)
        array_view[:] = values_arr

  def execute_impl(self, data_structure: nx.DataStructure, args: dict, message_handler: nx.IFilter.MessageHandler, should_cancel: nx.AtomicBoolProxy) -> nx.IFilter.ExecuteResult:
    cli_file_path: str = args[CliReaderFilter.CLI_FILE_PATH_KEY]
    output_edge_geom_path: nx.DataPath = args[CliReaderFilter.OUTPUT_EDGE_GEOM_PATH_KEY]
    output_edge_attrmat_name: str = args[CliReaderFilter.OUTPUT_EDGE_ATTRMAT_NAME_KEY]
    output_feature_attrmat_name: str = args[CliReaderFilter.OUTPUT_FEATURE_ATTRMAT_NAME_KEY]
    use_x_dimension_range: bool = args[CliReaderFilter.USE_X_DIMENSION_RANGE_KEY]
    use_y_dimension_range: bool = args[CliReaderFilter.USE_Y_DIMENSION_RANGE_KEY]
    use_z_dimension_range: bool = args[CliReaderFilter.USE_Z_DIMENSION_RANGE_KEY]
    out_of_bounds_behavior = CliReaderFilter.OutOfBoundsBehavior(args[CliReaderFilter.OUT_OF_BOUNDS_BEHAVIOR_KEY])
    min_max_x_coords: list = args[CliReaderFilter.MIN_MAX_X_COORDS_KEY]
    min_max_y_coords: list = args[CliReaderFilter.MIN_MAX_Y_COORDS_KEY]
    min_max_z_coords: list = args[CliReaderFilter.MIN_MAX_Z_COORDS_KEY]
    read_extra_metadata: bool = args[CliReaderFilter.READ_EXTRA_METADATA]
    create_edge_types_array: bool = args[CliReaderFilter.CREATE_EDGE_TYPE_ARRAY]

    bounding_box_coords = None
    if use_x_dimension_range or use_y_dimension_range or use_z_dimension_range:
      bounding_box_coords = [-sys.float_info.max, sys.float_info.max] * 3
    if use_x_dimension_range:
      bounding_box_coords[0:2] = min_max_x_coords
    if use_y_dimension_range:
      bounding_box_coords[2:4] = min_max_y_coords
    if use_z_dimension_range:
      bounding_box_coords[4:6] = min_max_z_coords

    try:
      result = self._parse_file(Path(cli_file_path),
                                out_of_bounds_behavior=out_of_bounds_behavior,
                                bounding_box=bounding_box_coords,
                                message_handler=message_handler)
      if result.invalid():
        return nx.Result(errors=result.errors)
      layer_features, layer_heights, hatch_labels = result.value
    except Exception as e:
      return nx.Result([nx.Error(-2010, f"An error occurred while parsing the CLI file '{cli_file_path}': {e}")])

    # Process the layer data to determine the start and end vertices, data arrays, and number of edges
    result = self._process_layers(layer_features, hatch_labels, message_handler)
    if result.invalid():
      return nx.Result(errors=result.errors)
    start_vertices, end_vertices, data_arrays, num_of_edges = result.value

    # Update the edge geometry
    self._update_edge_geometry(data_structure=data_structure,
                               start_vertices=start_vertices,
                               end_vertices=end_vertices,
                               data_arrays=data_arrays,
                               num_of_edges=num_of_edges,
                               output_edge_geom_path=output_edge_geom_path,
                               output_edge_attrmat_name=output_edge_attrmat_name,
                               output_feature_attrmat_name=output_feature_attrmat_name,
                               hatch_labels=hatch_labels,
                               read_extra_metadata=read_extra_metadata,
                               create_edge_types_array=create_edge_types_array,
                               message_handler=message_handler)

    # Filter is complete, return the results.
    return nx.Result()
  
  def _parse_file(self, full_path: Path, out_of_bounds_behavior: OutOfBoundsBehavior, bounding_box: list = None, message_handler: nx.IFilter.MessageHandler = None) -> Result:    
    layer_heights = None
    layer_features = None
    units = None
    
    #parse file lines
    with open(str(full_path), 'r') as file:
      while True:
        line = file.readline()
        if not line:
            break  # End of file
        
        line = line.strip()  # Remove surrounding whitespace
        if not line:
            continue  # Skip blank lines
        
        line = re.sub(r"//.*?//", "", line).strip()  # Remove comments
        if not line:
            continue  # Skip lines that become blank after removing comments
        
        if line.startswith("$$HEADERSTART"):
          if message_handler is not None:
            message_handler(nx.IFilter.Message(nx.IFilter.Message.Type.Info, 'Reading header data...'))
          units, hatch_labels = self._parse_header(file)
        if line.startswith("$$GEOMETRYSTART"):
          if message_handler is not None:
            message_handler(nx.IFilter.Message(nx.IFilter.Message.Type.Info, 'Reading geometry data...'))
          if units is None:
            return make_error_result(-8030, "Units have not been read and are needed when reading the geometry data.  Please make sure that a header with the $$HEADERSTART tag and containing the $$UNITS tag is read before this geometry.")
          result = self._parse_geometry(file, units, out_of_bounds_behavior, bounding_box, message_handler)
          if result.invalid():
            return Result(errors=result.errors)
          layer_features, layer_heights = result.value
    
    return Result(value=(layer_features, layer_heights, hatch_labels))
  
  def _parse_header(self, file):
    units = 1.  #default to 1, will be overwritten by anything read from the fileheader.
    hatch_labels = {}
    line = file.readline().strip()
    while not line.startswith("$$HEADEREND"):
      if line and line.startswith("$$UNITS"):
        key, val = line.split("/")
        units = float(val)
      elif line and line.startswith("$$LABEL"):
        key, val = line.split("/")
        hatch_id, hatch_desc = val.split(",")
        hatch_labels[int(hatch_id)] = hatch_desc.strip()
      line = file.readline()
        
    return units, hatch_labels

  def _parse_geometry_array_names(self, full_path: Path):
    array_names = []

    with open(str(full_path), 'r') as file:
      record = False
      num_of_labels = 0
      for line in file:
        line = re.sub(r"//.*?//", "", line).strip()  # Remove comments
        if not line:
          continue

        if line.startswith("$$LABEL"):
          num_of_labels = num_of_labels + 1

        if '$$GEOMETRYSTART' in line:
          record = True
          continue

        if line.startswith("$$HATCHES") or line.startswith("$HATCHES") or line.startswith("$$POLYLINE") or line.startswith("$POLYLINE"):
          break

        if '/' in line and record:
          tag, value_str = line.split("/", 1)
          tag = tag.replace('$', '')
          tag = tag.capitalize()
          if tag and tag not in array_names:
            array_names.append(tag)

    return array_names, num_of_labels

  def process_geometric_feature(self, feature_type, id, dir, n, coords, z_height, data, units, bounding_box, out_of_bounds_behavior, layer_counter, features):
    if bounding_box is not None:
        coords_result = self._filter_coordinates(coords, z_height, units, bounding_box, out_of_bounds_behavior)
        if coords_result.invalid():
            return Result(errors=coords_result.errors)
        coords = coords_result.value
        n = coords.size // 4
    start_xvals = coords[0::4] * units
    start_yvals = coords[1::4] * units
    end_xvals = coords[2::4] * units
    end_yvals = coords[3::4] * units
    new_feature = GeometricFeature(feature_type, layer_counter, z_height * units, data, id, n, start_xvals, start_yvals, end_xvals, end_yvals, dir)
    features.append(new_feature)

  def process_geometric_feature_line(self, read_line_func, feature_type, val, data, units, bounding_box, out_of_bounds_behavior, layer_counter, layer_heights, features) -> dict:
    id, dir, n, coords = read_line_func(val)
    z_height = layer_heights[-1]
    self.process_geometric_feature(feature_type, id, dir, n, coords, z_height, copy.copy(data), units, bounding_box, out_of_bounds_behavior, layer_counter, features)

  def process_extra_metadata_line(self, key: str, val: str, data: dict) -> dict:
    key = key.replace('$', '')
    key = key.capitalize()
    data[key] = val
    return data

  def _parse_geometry(self, file, units, out_of_bounds_behavior: OutOfBoundsBehavior, bounding_box: list = None, message_handler: nx.IFilter.MessageHandler = None) -> Result:
    layer_counter = -1 #initialize to -1, increment by one when finding the first layer
    layer_heights = []
    layer_features = []
    features = []
    data = {}
    
    last_message_time = 0
    
    #parse file lines
    line: str = file.readline()
    line = re.sub(r"//.*?//", "", line).strip()  # Remove comments
    while not line.startswith("$$GEOMETRYEND"):
      if not line:
        # Do nothing, read the next line
        line = file.readline()
        line = re.sub(r"//.*?//", "", line).strip()  # Remove comments
        continue

      key, val = line.split("/")
      if line.startswith("$$LAYER"):
        if layer_counter>=0:
          layer_features.append(features) #save the old stuff before starting the new stuff
        layer_counter += 1
        layer_heights.append(float(val))
        key = key.replace('$', '')
        key = key.capitalize()
        data[key] = val
        features = []  #reinitiate

        current_time = time.time()
        if message_handler is not None and (current_time - last_message_time) >= 1:
          message_handler(nx.IFilter.Message(nx.IFilter.Message.Type.Info, f"Number of layers found: {layer_counter + 1}"))
          last_message_time = current_time
      
      elif line.startswith("$$POLYLINE") or line.startswith("$POLYLINE"):
        self.process_geometric_feature_line(self._read_poly_line, 1, val, data, units, bounding_box, out_of_bounds_behavior, layer_counter, layer_heights, features)
      elif line.startswith("$$HATCHES") or line.startswith("$HATCHES"):
        self.process_geometric_feature_line(self._read_hatch_line, 0, val, data, units, bounding_box, out_of_bounds_behavior, layer_counter, layer_heights, features)
      else:
        data = self.process_extra_metadata_line(key, val, data)

      line = file.readline()
      line = re.sub(r"//.*?//", "", line).strip()  # Remove comments
    
    message_handler(nx.IFilter.Message(nx.IFilter.Message.Type.Info, f"Number of layers found: {layer_counter + 1}"))

    #convert heights into an array, and scale by the units value
    layer_heights = units * np.array(layer_heights)

    # Save the last layer
    if features:
      layer_features.append(features)

    return Result(value=(layer_features, layer_heights))

  def _read_poly_line(self, line: str) -> Tuple[int, int, int, np.ndarray]:
    vals = line.split(",")

    # Set up coord_array so that it contains start & end coordinate pairs, not just a list of coordinates.
    # This will make the rest of the code work for polylines like it does for hatches.
    coord_array = np.array(list(map(float, vals[3:])))
    reshaped_arr = coord_array.reshape(int(len(coord_array) / 2), 2)
    repeated_arr = np.repeat(reshaped_arr[1:-1], 2, axis=0)
    coord_array = np.concatenate(([reshaped_arr[0]], repeated_arr, [reshaped_arr[-1]]))
    coord_array = coord_array.ravel()

    return int(vals[0]), int(vals[1]), int(len(coord_array) / 4), coord_array

  def _read_hatch_line(self, line: str) -> Tuple[int, int, np.ndarray]:
    vals = line.split(",")
    return int(vals[0]), None, int(vals[1]), np.array(list(map(float, vals[2:])))

  def _filter_coordinates(self, coords: np.ndarray, z_height: float, units: int, bounding_box: list, out_of_bounds_behavior: OutOfBoundsBehavior) -> Result[np.array]:
    x_min, x_max, y_min, y_max, z_min, z_max = bounding_box
    
    coords_output = coords.copy()
    coords_output = coords_output.reshape((coords.size // 4, 4))

    # Extracting x and y coordinates from the numpy array
    x1 = coords_output[:, 0] * units
    y1 = coords_output[:, 1] * units
    x2 = coords_output[:, 2] * units
    y2 = coords_output[:, 3] * units
    z = z_height * units

    # Check for violations of bounding box condition
    x1_inside = (x1 >= x_min) & (x1 <= x_max)
    y1_inside = (y1 >= y_min) & (y1 <= y_max)
    x2_inside = (x2 >= x_min) & (x2 <= x_max)
    y2_inside = (y2 >= y_min) & (y2 <= y_max)
    z_inside = (z >= z_min) & (z <= z_max)

    if out_of_bounds_behavior == CliReaderFilter.OutOfBoundsBehavior.FilterError:
      # Check for cases where one vertex is inside the bounding box while the other is outside the bounding box
      if np.any(x1_inside & y1_inside & z_inside & (~x2_inside | ~y2_inside)):
        idx = np.where(x1_inside & y1_inside & z_inside & (~x2_inside | ~y2_inside))[0][0]
        return make_error_result(code=-1100, message=f"Point ({x1[idx]}, {y1[idx]}, {z}) is inside the bounding box while point ({x2[idx]}, {y2[idx]}, {z}) is outside.")
      if np.any(x2_inside & y2_inside & z_inside & (~x1_inside | ~y1_inside)):
        idx = np.where(x2_inside & y2_inside & z_inside & (~x1_inside | ~y1_inside))[0][0]
        return make_error_result(code=-1101, message=f"Point ({x2[idx]}, {y2[idx]}, {z}) is inside the bounding box while point ({x1[idx]}, {y1[idx]}, {z}) is outside.")
    
    if out_of_bounds_behavior == CliReaderFilter.OutOfBoundsBehavior.InterpolateOutsideVertex:
      for i in range(coords_output.shape[0]):
        if (x1_inside[i] & y1_inside[i] & z_inside) and not (x2_inside[i] & y2_inside[i]):
          x2[i], y2[i] = self._interpolate_outside_vertex((x1[i], y1[i]), (x2[i], y2[i]), x_min, x_max, y_min, y_max)
          x2_inside[i] = True
          y2_inside[i] = True
        elif (x2_inside[i] & y2_inside[i] & z_inside) and not (x1_inside[i] & y1_inside[i]):
          x1[i], y1[i] = self._interpolate_outside_vertex((x2[i], y2[i]), (x1[i], y1[i]), x_min, x_max, y_min, y_max)
          x1_inside[i] = True
          y1_inside[i] = True

      coords_output[:, 0] = x1 / units
      coords_output[:, 1] = y1 / units
      coords_output[:, 2] = x2 / units
      coords_output[:, 3] = y2 / units

    # Creating boolean masks for points inside the bounding box
    mask = x1_inside & y1_inside & x2_inside & y2_inside & z_inside

    # Filtering numpy array based on the mask
    coords_output = coords_output[mask]

    coords_output = coords_output.reshape((np.prod(coords_output.shape)))
    return Result(value=coords_output)
  
  def _interpolate_outside_vertex(self, inside_vertex, outside_vertex, x_min, x_max, y_min, y_max):
    """
    Adjusts the outside vertex along a line segment so that it remains within the specified bounding box.

    Parameters:
    - inside_vertex (tuple): The inside vertex (x1, y1) of the line segment.
    - outside_vertex (tuple): The outside vertex (x2, y2) of the line segment that is outside the bounding box.
    - x_min (float): Minimum x-coordinate of the bounding box.
    - x_max (float): Maximum x-coordinate of the bounding box.
    - y_min (float): Minimum y-coordinate of the bounding box.
    - y_max (float): Maximum y-coordinate of the bounding box.

    Returns:
    - tuple: Adjusted coordinates of outside_vertex such that the entire line segment lies within the bounding box.
    """

    x1, y1 = inside_vertex
    x2, y2 = outside_vertex

    def clip(t_min, t_max, v1, v2, vmin, vmax):
      """
      Adjusts the range of t along a line segment to ensure the segment remains within bounds.

      Parameters:
      - t_min (float): Current minimum of parameter t for the line segment.
      - t_max (float): Current maximum of parameter t for the line segment.
      - v1 (float): Coordinate (x or y) of the start point of the line segment.
      - v2 (float): Coordinate (x or y) of the end point of the line segment.
      - vmin (float): Minimum allowable value for the coordinate.
      - vmax (float): Maximum allowable value for the coordinate.

      Returns:
      - (float, float): Updated minimum and maximum values of t after applying the clipping based on bounds.
      """
      if v1 < v2:  # Moving towards increasing coordinate
        if v2 > vmax:
          t_max = min(t_max, (vmax - v1) / (v2 - v1))
        if v2 < vmin:
          t_min = max(t_min, (vmin - v1) / (v2 - v1))
      else:  # Moving towards decreasing coordinate
        if v2 < vmin:
          t_max = min(t_max, (vmin - v1) / (v2 - v1))
        if v2 > vmax:
          t_min = max(t_min, (vmax - v1) / (v2 - v1))
      return t_min, t_max

    t_min, t_max = 0.0, 1.0

    t_min, t_max = clip(t_min, t_max, x1, x2, x_min, x_max)
    t_min, t_max = clip(t_min, t_max, y1, y2, y_min, y_max)

    if t_min > t_max:
        return x2, y2  # No intersection, return original outside vertex

    t = t_max  # Use t_max to get the intersection point

    new_x = x1 + t * (x2 - x1)
    new_y = y1 + t * (y2 - y1)

    return new_x, new_y
