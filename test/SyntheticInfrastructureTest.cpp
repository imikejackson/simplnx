//
// Created by Michael Jackson on 3/12/25.
//

#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/DataStructure/IStatsData.hpp"
#include "simplnx/unit_test/simplnx_test_dirs.hpp"

#include <catch2/catch.hpp>

#include <string>

using namespace nx::core;

TEST_CASE("Simplnx::SyntheticTests", "[Simplnx][SyntheticTests]")
{

  DataStructure dataStructure;
  auto statsDataPtr = StatsDataArray::Create(dataStructure, "test stats data array");

  auto primaryStatsData = std::make_unique<PrimaryStatsData>();

  statsDataPtr->addShape(std::move(primaryStatsData));
}
