// Copyright 2020 The XLS Authors
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "xls/noc/config/router_config_proto_builder.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace xls::noc {

// Test the field values of a router.
TEST(RouterConfigBuilderTest, FieldValues) {
  const char* kName = "Test";
  const char* kRouterPortName = "RouterPort";
  RouterConfigProto proto;
  RouterConfigProtoBuilder builder(&proto);
  builder.WithName(kName);
  PortConfigProtoBuilder port_config_builder =
      builder.WithPort(kRouterPortName);
  port_config_builder.AsOutputDirection();
  RoutingSchemeConfigProtoBuilder routing_scheme_config_builder =
      builder.GetRoutingSchemeConfigProtoBuilder();
  routing_scheme_config_builder.WithDistributedRoutingEntry(
      {"NetworkPort", "VC0"}, {"RouterPort", "VC0"});
  ArbiterSchemeConfigProtoBuilder arbiter_scheme_config_builder =
      builder.GetArbiterSchemeConfigProtoBuilder();
  std::vector<PortVirtualChannelTuple> priority_list;
  priority_list.push_back({"RouterInputPort", "VC0"});
  arbiter_scheme_config_builder.WithPriorityEntry(kRouterPortName,
                                                  priority_list);

  EXPECT_TRUE(proto.has_name());
  EXPECT_EQ(proto.ports_size(), 1);
  EXPECT_THAT(proto.ports(0).name(), kRouterPortName);
  EXPECT_TRUE(proto.has_routing_scheme());
  EXPECT_TRUE(proto.has_arbiter_scheme());
  EXPECT_THAT(proto.name(), kName);
}

}  // namespace xls::noc
