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

test tuple_index {
  // Perform tuple indexing without explicitly annotating the type of the index.
  let _ = assert_eq((u32:42, u4:7)[1], u4:7);
  let _ = assert_eq((u32:42, u4:7)[0], u32:42);
  ()
}

// TODO(leary): 2020-11-09 Sample that will only work when we add unifying type
// inference (the type variable that results from the indexing operation must
// unify with the uN[R] result).
//
//fn index_parametric<N: u32, R: u32>(t: (u32, u4)) -> uN[R] {
//  t[N]
//}
//
//test tuple_index_parametric {
//  let t = (u32:42, u4:7);
//  let _ = assert_eq(index_parametric<u32:0, u32:32>(t), u32:42);
//  let _ = assert_eq(index_parametric<u32:1, u32:4>(t), u4:7);
//  ()
//}
