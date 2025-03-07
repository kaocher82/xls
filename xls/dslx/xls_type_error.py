# Lint as: python3
#
# Copyright 2020 The XLS Authors
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Defines the type error that DSLX can produce during type checking."""

from typing import Text, Optional, Tuple

from xls.dslx.python import cpp_pos
from xls.dslx.python.cpp_concrete_type import ConcreteType
from xls.dslx.span import PositionalError


class ArgCountMismatchError(PositionalError):
  """Raised when argument count != parameter count in an invocation."""

  def __init__(self, span: cpp_pos.Span, arg_types: Tuple[ConcreteType, ...],
               param_count: int, param_types: Optional[Tuple[ConcreteType,
                                                             ...]],
               suffix: Optional[Text]):
    self.arg_types = arg_types
    self.param_types = param_types
    parameters_str = ' parameters: [{}];'.format(', '.join(
        str(p) for p in param_types)) if param_types else ''
    msg = ('Expected {} parameter(s) but got {} argument(s);{} '
           'arguments: [{}]').format(param_count, len(arg_types),
                                     parameters_str,
                                     ', '.join(str(a) for a in arg_types))
    super(ArgCountMismatchError, self).__init__(msg, span)
