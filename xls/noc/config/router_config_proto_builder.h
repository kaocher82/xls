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

#ifndef XLS_NOC_CONFIG_ROUTER_CONFIG_PROTO_BUILDER_H_
#define XLS_NOC_CONFIG_ROUTER_CONFIG_PROTO_BUILDER_H_

#include "xls/noc/config/arbiter_scheme_config_proto_builder.h"
#include "xls/noc/config/network_config.pb.h"
#include "xls/noc/config/port_config_proto_builder.h"
#include "xls/noc/config/routing_scheme_config_proto_builder.h"

namespace xls::noc {

// A builder for constructing a router configuration proto.
class RouterConfigProtoBuilder {
 public:
  // proto cannot be nullptr.
  explicit RouterConfigProtoBuilder(RouterConfigProto* proto)
      : proto_(XLS_DIE_IF_NULL(proto)) {}

  // Adds the name of the router.
  RouterConfigProtoBuilder& WithName(absl::string_view name);

  // Adds a port to this object. Returns its builder.
  PortConfigProtoBuilder WithPort(absl::string_view name);

  // Returns the routing scheme configuration builder of this builder.
  RoutingSchemeConfigProtoBuilder GetRoutingSchemeConfigProtoBuilder();

  // Returns the arbiter scheme configuration builder of this builder.
  ArbiterSchemeConfigProtoBuilder GetArbiterSchemeConfigProtoBuilder();

 private:
  RouterConfigProto* proto_;
};

}  // namespace xls::noc

#endif  // XLS_NOC_CONFIG_ROUTER_CONFIG_PROTO_BUILDER_H_
