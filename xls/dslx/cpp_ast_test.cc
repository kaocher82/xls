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

#include "xls/dslx/cpp_ast.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "xls/common/status/matchers.h"

namespace xls::dslx {
namespace {

using status_testing::IsOkAndHolds;
using status_testing::StatusIs;
using testing::HasSubstr;

TEST(CppAst, ModuleWithConstant) {
  Module m("test");
  const Span fake_span;
  Number* number = m.Make<Number>(fake_span, std::string("42"),
                                  NumberKind::kOther, /*type=*/nullptr);
  NameDef* name_def = m.Make<NameDef>(fake_span, std::string("MOL"), nullptr);
  ConstantDef* constant_def =
      m.Make<ConstantDef>(fake_span, name_def, number, /*is_public=*/false);
  name_def->set_definer(constant_def);
  m.AddTop(constant_def);

  EXPECT_EQ(m.ToString(), "const MOL = 42;");
}

TEST(CppAst, GetNumberAsInt64) {
  struct Example {
    std::string text;
    int64 want;
  } kCases[] = {
      {"0b0", 0},
      {"0b1", 1},
      {"0b10", 2},
      {"0b11", 3},
      {"0b100", 4},
      {"0b1000", 8},
      {"0b1011", 11},
      {"0b1_1000", 24},
      {"0b1_1001", 25},
      {"0b1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_1111_"
       "1111_1111_1111",
       -1},
  };
  Module m("test");
  auto make_num = [&m](std::string text) {
    const Span fake_span;
    return m.Make<Number>(fake_span, text, NumberKind::kOther,
                          /*type=*/nullptr);
  };
  for (const Example& example : kCases) {
    EXPECT_THAT(make_num(example.text)->GetAsInt64(),
                IsOkAndHolds(example.want));
  }

  EXPECT_THAT(make_num("0b")->GetAsInt64(),
              StatusIs(absl::StatusCode::kInvalidArgument,
                       HasSubstr("Could not convert 0b to a number")));
}

}  // namespace
}  // namespace xls::dslx
