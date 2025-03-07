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

#ifndef XLS_DSLX_INTERP_VALUE_H_
#define XLS_DSLX_INTERP_VALUE_H_

#include "xls/dslx/cpp_ast.h"
#include "xls/ir/bits.h"

namespace xls::dslx {

#define XLS_DSLX_BUILTIN_EACH(X)     \
  X("add_with_carry", kAddWithCarry) \
  X("and_reduce", kAndReduce)        \
  X("assert_eq", kAssertEq)          \
  X("assert_lt", kAssertLt)          \
  X("bit_slice", kBitSlice)          \
  X("clz", kClz)                     \
  X("ctz", kCtz)                     \
  X("enumerate", kEnumerate)         \
  X("fail!", kFail)                  \
  X("map", kMap)                     \
  X("one_hot", kOneHot)              \
  X("one_hot_sel", kOneHotSel)       \
  X("or_reduce", kOrReduce)          \
  X("range", kRange)                 \
  X("rev", kRev)                     \
  X("select", kSelect)               \
  X("sge", kSGe)                     \
  X("sgt", kSGt)                     \
  X("signex", kSignex)               \
  X("sle", kSLe)                     \
  X("slice", kSlice)                 \
  X("slt", kSLt)                     \
  X("trace", kTrace)                 \
  X("update", kUpdate)               \
  X("xor_reduce", kXorReduce)

// Enum that represents all the DSLX builtin functions.
//
// Functions can be held in values, either as user defined ones or builtin ones
// (represented via this enumerated value).
enum class Builtin {
#define ENUMIFY(__str, __enum, ...) __enum,
  XLS_DSLX_BUILTIN_EACH(ENUMIFY)
#undef ENUMIFY
};

absl::StatusOr<Builtin> BuiltinFromString(absl::string_view name);

std::string BuiltinToString(Builtin builtin);

static const Builtin kAllBuiltins[] = {
#define ELEMIFY(__str, __enum, ...) Builtin::__enum,
    XLS_DSLX_BUILTIN_EACH(ELEMIFY)
#undef ELEMIFY
};

// Tags a value to denote its payload.
//
// Note this goes beyond InterpValue::Payload annotating things like whether the
// bits should be interpreted as signed or unsigned, which can change the
// behavior of interpreted operators like '<'.
enum class InterpValueTag {
  kUBits,
  kSBits,
  kTuple,
  kArray,
  kEnum,
  kFunction,
};

std::string TagToString(InterpValueTag tag);

inline std::ostream& operator<<(std::ostream& os, InterpValueTag tag) {
  os << TagToString(tag);
  return os;
}

// A DSLX interpreter value (variant), with InterpValueTag as a discriminator.
class InterpValue {
 public:
  struct UserFnData {
    std::shared_ptr<Module> module;
    Function* function;
  };
  using FnData = absl::variant<Builtin, UserFnData>;

  // Factories

  // TODO(leary): 2020-10-18 Port to be consistent with xls/ir/bits.h
  static InterpValue MakeUBits(int64 bit_count, int64 value);
  static InterpValue MakeSBits(int64 bit_count, int64 value);

  static InterpValue MakeNil() { return MakeTuple({}); }
  static InterpValue MakeU32(uint32 value) {
    return MakeUBits(/*bit_count=*/32, value);
  }
  static InterpValue MakeU64(uint64 value) {
    return MakeUBits(/*bit_count=*/64, value);
  }
  static InterpValue MakeEnum(const Bits& bits, EnumDef* type) {
    return InterpValue(InterpValueTag::kEnum, bits, type);
  }
  static InterpValue MakeTuple(std::vector<InterpValue> members);
  static absl::StatusOr<InterpValue> MakeArray(
      std::vector<InterpValue> elements);
  static InterpValue MakeBool(bool value) { return MakeUBits(1, value); }

  static absl::StatusOr<InterpValue> MakeBits(InterpValueTag tag, Bits bits) {
    if (tag != InterpValueTag::kUBits && tag != InterpValueTag::kSBits) {
      return absl::InvalidArgumentError(
          "Bits tag is required to make a bits Value; got: " +
          TagToString(tag));
    }
    return InterpValue(tag, std::move(bits));
  }

  static InterpValue MakeFunction(Builtin b) {
    return InterpValue(InterpValueTag::kFunction, b);
  }
  static InterpValue MakeFunction(FnData fn_data) {
    return InterpValue(InterpValueTag::kFunction, std::move(fn_data));
  }

  // Queries

  bool IsTuple() const { return tag_ == InterpValueTag::kTuple; }
  bool IsArray() const { return tag_ == InterpValueTag::kArray; }
  bool IsNilTuple() const { return IsTuple() && GetLength().value() == 0; }
  bool IsUBits() const { return tag_ == InterpValueTag::kUBits; }
  bool IsSBits() const { return tag_ == InterpValueTag::kSBits; }
  bool IsBits() const { return IsUBits() || IsSBits(); }
  bool IsEnum() const { return tag_ == InterpValueTag::kEnum; }
  bool IsFunction() const { return tag_ == InterpValueTag::kFunction; }
  bool IsBuiltinFunction() const {
    return IsFunction() && absl::holds_alternative<Builtin>(GetFunctionOrDie());
  }

  bool IsTraceBuiltin() const {
    return IsBuiltinFunction() &&
           absl::get<Builtin>(GetFunctionOrDie()) == Builtin::kTrace;
  }

  bool IsFalse() const { return IsBool() && GetBitsOrDie().IsZero(); }
  bool IsTrue() const { return IsBool() && GetBitsOrDie().IsAllOnes(); }

  bool IsSigned() const {
    XLS_CHECK(IsBits() || IsEnum());
    return IsSBits() || (IsEnum() && type()->signedness().value());
  }

  // Note that this equality for bit value holding tags ignores the tag and just
  // checks bit pattern equivalence. We expect that the type checker will
  // validate the comparisons are done properly between types before the
  // interpreter performs these operations.
  bool Eq(const InterpValue& other) const;
  bool Ne(const InterpValue& other) const { return !Eq(other); }

  // Various operations -- tag is not applicable for the operation, a status
  // error should be returned.

  absl::StatusOr<int64> GetBitValueCheckSign() const;

  absl::StatusOr<int64> GetBitCount() const;
  absl::StatusOr<uint64> GetBitValueUint64() const;
  absl::StatusOr<int64> GetBitValueInt64() const;

  absl::StatusOr<int64> GetLength() const {
    if (IsTuple() || IsArray()) {
      return GetValuesOrDie().size();
    }
    return absl::InvalidArgumentError("Invalid tag for length query: " +
                                      TagToString(tag_));
  }

  absl::StatusOr<InterpValue> ZeroExt(int64 new_bit_count) const;
  absl::StatusOr<InterpValue> SignExt(int64 new_bit_count) const;
  absl::StatusOr<InterpValue> Concat(const InterpValue& other) const;
  absl::StatusOr<InterpValue> AddWithCarry(const InterpValue& other) const;

  absl::StatusOr<InterpValue> Add(const InterpValue& other) const;
  absl::StatusOr<InterpValue> Sub(const InterpValue& other) const;

  absl::StatusOr<InterpValue> Mul(const InterpValue& other) const;
  absl::StatusOr<InterpValue> Shll(const InterpValue& other) const;
  absl::StatusOr<InterpValue> Shrl(const InterpValue& other) const;
  absl::StatusOr<InterpValue> Shra(const InterpValue& other) const;
  absl::StatusOr<InterpValue> BitwiseXor(const InterpValue& other) const;
  absl::StatusOr<InterpValue> BitwiseNegate() const;
  absl::StatusOr<InterpValue> BitwiseOr(const InterpValue& other) const;
  absl::StatusOr<InterpValue> BitwiseAnd(const InterpValue& other) const;
  absl::StatusOr<InterpValue> ArithmeticNegate() const;
  absl::StatusOr<InterpValue> FloorDiv(const InterpValue& other) const;
  absl::StatusOr<InterpValue> Index(const InterpValue& other) const;
  absl::StatusOr<InterpValue> Update(const InterpValue& index,
                                     const InterpValue& value) const;
  absl::StatusOr<InterpValue> Slice(const InterpValue& start,
                                    const InterpValue& length) const;
  absl::StatusOr<InterpValue> Flatten() const;
  absl::StatusOr<InterpValue> OneHot(bool lsb_prio) const;

  absl::StatusOr<InterpValue> Lt(const InterpValue& other) const;
  absl::StatusOr<InterpValue> Le(const InterpValue& other) const;
  absl::StatusOr<InterpValue> Gt(const InterpValue& other) const;
  absl::StatusOr<InterpValue> Ge(const InterpValue& other) const;

  // Performs the signed comparison defined by "method".
  //
  // "method" should be a value in the set {slt, sle, sgt, sge} or an
  // InvalidArgumentError is returned.
  absl::StatusOr<InterpValue> SCmp(const InterpValue& other,
                                   absl::string_view method);

  // "humanize" indicates whether bits-based values /within containers/ should
  // not have their leading widths printed -- this can obscure readability.
  std::string ToString(bool humanize = false) const;
  std::string ToHumanString() const { return ToString(/*humanize=*/true); }

  InterpValueTag tag() const { return tag_; }

  absl::StatusOr<const std::vector<InterpValue>*> GetValues() const {
    if (!absl::holds_alternative<std::vector<InterpValue>>(payload_)) {
      return absl::InvalidArgumentError("Value does not hold element values");
    }
    return &absl::get<std::vector<InterpValue>>(payload_);
  }
  const std::vector<InterpValue>& GetValuesOrDie() const {
    return absl::get<std::vector<InterpValue>>(payload_);
  }
  absl::StatusOr<const FnData*> GetFunction() const {
    if (!absl::holds_alternative<FnData>(payload_)) {
      return absl::InvalidArgumentError("Value does not hold function data");
    }
    return &absl::get<FnData>(payload_);
  }
  const FnData& GetFunctionOrDie() const { return *GetFunction().value(); }
  absl::StatusOr<Bits> GetBits() const;
  const Bits& GetBitsOrDie() const { return absl::get<Bits>(payload_); }

  // For enum values, the enum that the bit pattern is interpreted by is
  // referred to by the interpreter value.
  EnumDef* type() const { return type_; }

  // Note: different from IsBits() which is checking whether the tag is sbits or
  // ubits; this is checking whether there are bits in the payload, which would
  // apply to enum values as well.
  bool HasBits() const { return absl::holds_alternative<Bits>(payload_); }
  bool HasValues() const {
    return absl::holds_alternative<std::vector<InterpValue>>(payload_);
  }

 private:
  friend struct InterpValuePickler;

  using Payload = absl::variant<Bits, std::vector<InterpValue>, FnData>;

  InterpValue(InterpValueTag tag, Payload payload, EnumDef* type = nullptr)
      : tag_(tag), payload_(std::move(payload)), type_(type) {}

  using CompareF = bool (*)(const Bits& lhs, const Bits& rhs);

  bool IsBool() const { return IsBits() && GetBitCount().value() == 1; }

  // Helper for various comparisons.
  static absl::StatusOr<InterpValue> Compare(const InterpValue& lhs,
                                             const InterpValue& rhs,
                                             CompareF ucmp, CompareF scmp);

  InterpValueTag tag_;
  Payload payload_;

  // For enum values (where we refer to the enum definition) we keep the module
  // holding the enum AST node alive via the shared_ptr. When the DSL is fully
  // migrated to C++ we can use lifetime assumptions as is more typical in C++
  // land, this is done for Python interop.
  EnumDef* type_ = nullptr;
};

}  // namespace xls::dslx

#endif  // XLS_DSLX_INTERP_VALUE_H_
