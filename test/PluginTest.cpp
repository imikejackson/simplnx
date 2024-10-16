#include "simplnx/Core/Application.hpp"
#include "simplnx/DataStructure/DataStructure.hpp"
#include "simplnx/Filter/FilterHandle.hpp"
#include "simplnx/Filter/IFilter.hpp"
#include "simplnx/Utilities/StringUtilities.hpp"
#include "simplnx/unit_test/simplnx_test_dirs.hpp"

#include <catch2/catch.hpp>

#include <sstream>
#include <string>

using namespace nx::core;

namespace
{
constexpr Uuid k_TestOnePluginId = *Uuid::FromString("01ff618b-781f-4ac0-b9ac-43f26ce1854f");
constexpr Uuid k_TestFilterId = *Uuid::FromString("5502c3f7-37a8-4a86-b003-1c856be02491");
const FilterHandle k_TestFilterHandle(k_TestFilterId, k_TestOnePluginId);

constexpr Uuid k_TestTwoPluginId = *Uuid::FromString("05cc618b-781f-4ac0-b9ac-43f33ce1854e");
constexpr Uuid k_Test2FilterId = *Uuid::FromString("ad9cf22b-bc5e-41d6-b02e-bb49ffd12c04");
const FilterHandle k_Test2FilterHandle(k_Test2FilterId, k_TestTwoPluginId);
} // namespace

TEST_CASE("Test Loading Plugins")
{
  auto app = Application::GetOrCreateInstance();
  app->loadPlugins(unit_test::k_BuildDir.view());

  auto* filterListPtr = Application::Instance()->getFilterList();
  const auto& filterHandles = filterListPtr->getFilterHandles();
  auto plugins = filterListPtr->getLoadedPlugins();

  if(plugins.size() != SIMPLNX_PLUGIN_COUNT)
  {
    std::cout << "Incorrect number of plugins were loaded.\n"
              << "Expected: " << SIMPLNX_PLUGIN_COUNT << "\nLoaded: " << plugins.size() << "\nLoaded Plugins are:\n";
    for(auto const& plugin : plugins)
    {
      std::cout << plugin->getName() << "\n";
    }
  }

  REQUIRE(plugins.size() == SIMPLNX_PLUGIN_COUNT);
  REQUIRE(filterHandles.size() >= 2);

  DataStructure dataStructure;
  {

    IFilter::UniquePointer filter = filterListPtr->createFilter(k_TestFilterHandle);
    REQUIRE(filter != nullptr);
    REQUIRE(filter->humanName() == "Test Filter");
    filter->execute(dataStructure, {});
  }
  {
    IFilter::UniquePointer filter2 = filterListPtr->createFilter(k_Test2FilterHandle);
    REQUIRE(filter2 != nullptr);
    REQUIRE(filter2->humanName() == "Test Filter 2");
    filter2->execute(dataStructure, {});
  }
}

TEST_CASE("Test Singleton")
{
  auto app = Application::GetOrCreateInstance();
  app->loadPlugins(unit_test::k_BuildDir.view());

  REQUIRE(app != nullptr);

  auto* filterListPtr = app->getFilterList();
  const auto& filterHandles = filterListPtr->getFilterHandles();
  auto plugins = filterListPtr->getLoadedPlugins();

  // Check plugins were loaded
  if(plugins.size() != SIMPLNX_PLUGIN_COUNT)
  {
    std::cout << "Incorrect number of plugins were loaded.\n"
              << "Expected: " << SIMPLNX_PLUGIN_COUNT << "\nLoaded: " << plugins.size() << "\nLoaded Plugins are:\n";
    for(auto const& plugin : plugins)
    {
      std::cout << plugin->getName() << "\n";
    }
  }
  REQUIRE(plugins.size() == SIMPLNX_PLUGIN_COUNT);

  // Check filters loaded
  REQUIRE(filterHandles.size() >= 2);

  // Create and execute filters
  DataStructure dataStructure;
  {
    IFilter::UniquePointer filter = filterListPtr->createFilter(k_TestFilterHandle);
    REQUIRE(filter != nullptr);
    REQUIRE(filter->humanName() == "Test Filter");
    filter->execute(dataStructure, {});
  }
  {
    IFilter::UniquePointer filter2 = filterListPtr->createFilter(k_Test2FilterHandle);
    REQUIRE(filter2 != nullptr);
    REQUIRE(filter2->humanName() == "Test Filter 2");
    filter2->execute(dataStructure, {});
  }

  Application::DeleteInstance();
  REQUIRE(Application::Instance() == nullptr);
}
