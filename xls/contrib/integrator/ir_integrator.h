// Copyright 2020 Google LLC
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

#ifndef XLS_INTEGRATOR_IR_INTEGRATOR_H_
#define XLS_INTEGRATOR_IR_INTEGRATOR_H_

#include "absl/status/statusor.h"
#include "xls/ir/function.h"
#include "xls/ir/package.h"

namespace xls {

class IntegrationOptions {
 public:
  // Used to specify different integration algorithms.
  enum class Algorithm {
    kBasicIntegrationAlgorithm,
  };

  // Which algorithm to use to merge functions.
  IntegrationOptions& algorithm(Algorithm value) {
    algorithm_ = value;
    return *this;
  }
  Algorithm algorithm() const { return algorithm_; }

  // Whether we can program individual muxes with unique select signals
  // or if we can configure the entire graph to match one of the input
  // functions using a single select signal.
  IntegrationOptions& unique_select_signal_per_mux(bool value) {
    unique_select_signal_per_mux_ = value;
    return *this;
  }
  bool unique_select_signal_per_mux() const {
    return unique_select_signal_per_mux_;
  }

 private:
  bool unique_select_signal_per_mux_ = false;
  Algorithm algorithm_ = Algorithm::kBasicIntegrationAlgorithm;
};

std::ostream& operator<<(std::ostream& os,
                         const IntegrationOptions::Algorithm& alg) {
  switch (alg) {
    case IntegrationOptions::Algorithm::kBasicIntegrationAlgorithm:
      os << "kBasicIntegrationAlgorithm";
      break;
  }
  return os;
}

// Class that represents an integration function i.e. a function combining the
// IR of other functions. This class tracks which original function nodes are
// mapped to which integration function nodes. It also provides some utilities
// that are useful for constructing the integrated function.
// TODO(jbaileyhandle): Consider breaking this up into 2 or more smaller
// classes (e.g. could have a class for handling muxes, another for tracking
// mappings etc). Should move these classes and their tests into a new file(s).
class IntegrationFunction {
 public:
  IntegrationFunction(const IntegrationFunction& other) = delete;
  void operator=(const IntegrationFunction& other) = delete;

  // Create an IntegrationFunction object that is empty expect for
  // parameters. Each initial parameter of the function is a tuple
  // which holds inputs corresponding to the parameters of one
  // of the source_functions.
  static absl::StatusOr<std::unique_ptr<IntegrationFunction>>
  MakeIntegrationFunctionWithParamTuples(
      Package* package, absl::Span<const Function* const> source_functions,
      const IntegrationOptions& options = IntegrationOptions(),
      std::string function_name = "IntegrationFunction");

  // Create a tuple of the nodes that are the map targets of
  // the return values of the source functions. Set this value
  // as the return value of the integration function and return the tuple.
  // Should not be called if a return value is already set for
  // the integration function.
  absl::StatusOr<Node*> MakeTupleReturnValue();

  Function* function() const { return function_.get(); }
  Node* global_mux_select() const { return global_mux_select_param_; }

  // Estimate the cost of inserting the node 'to_insert'.
  absl::StatusOr<float> GetInsertNodeCost(const Node* to_insert);

  // Add the external node 'to_insert' into the integration function. Mapped
  // operands are automatically discovered and connected.
  absl::StatusOr<Node*> InsertNode(Node* to_insert);

  // Estimate the cost of merging node_a and node_b. If nodes cannot be
  // merged, not value is returned.
  absl::StatusOr<std::optional<int64>> GetMergeNodesCost(const Node* node_a,
                                                         const Node* node_b);

  // Merge node_a and node_b. Operands are automatically multiplexed.
  // Returns the nodes that node_a and node_b map to after merging (vector
  // contains a single node if they map to the same node).
  absl::StatusOr<std::vector<Node*>> MergeNodes(Node* node_a, Node* node_b);

  enum class UnificationChange {
    kNoChange,
    kNewMuxAdded,
    kExistingMuxCasesModified,
  };
  struct UnifiedNode {
    Node* node;
    UnificationChange change;
  };
  // For the integration function nodes node_a and node_b,
  // returns a UnifiedNode. UnifiedNode.node points to a single integration
  // function node that combines the two nodes. This may involve adding a mux
  // and a parameter to serve as the mux select signal.
  // UnifiedNode.change will indicate if the ir graph is unchaged, if a new
  // mux was added by this call, or if an existing mux was modified by this
  // call.
  absl::StatusOr<UnifiedNode> UnifyIntegrationNodes(Node* node_a, Node* node_b);

  // Return a UnifiedOperands struct in which the 'operands' vector holds nodes
  // where each node unifies the corresponding operands of 'node_a' and
  // 'node_b'.  The 'changed_muxes' field lists all muxes created by or modified
  // by this call.
  struct UnifiedOperands {
    std::vector<Node*> operands;
    std::vector<UnifiedNode> changed_muxes;
  };
  absl::StatusOr<UnifiedOperands> UnifyNodeOperands(const Node* node_a,
                                                    const Node* node_b);

  // Return the case indexes which are in active use for a mux
  // whose selector is global_mux_select_param_;
  absl::StatusOr<const std::set<int64>*> GetGlobalMuxOccupiedCaseIndexes(
      const Node* node) const;

  // Return the case indices which were most recently added to a mux
  // whose selector is global_mux_select_param_;
  absl::StatusOr<const std::set<int64>*> GetGlobalMuxLastCaseIndexesAdded(
      const Node* node) const;

  // Returns how many muxes whose selector is global_mux_select_param_
  // we track metadata for.
  int64 GetNumberOfGlobalMuxesTracked() const;

  // For a mux produced by UnifyIntegrationNodes, undoes the previous
  // call to UnifyIntegrationNodes. Also updates
  // internal unification / mux book-keeping accordingly. If the mux was
  // modified / replaced, the new mux is returned. If the mux was removed,
  // nullptr is returned.
  absl::StatusOr<Node*> DeUnifyIntegrationNodes(Node* node);

  // Declares that node 'source' from a source function maps
  // to node 'map_target' in the integrated_function.
  absl::Status SetNodeMapping(const Node* source, Node* map_target);

  // Returns the integrated node that 'original' maps to, if it
  // exists. Otherwise, return an error status.
  absl::StatusOr<Node*> GetNodeMapping(const Node* original) const;

  // Returns the original nodes that map to 'map_target' in the integrated
  // function.
  absl::StatusOr<const absl::flat_hash_set<const Node*>*> GetNodesMappedToNode(
      const Node* map_target) const;

  // Returns the source function index for the given node.
  absl::StatusOr<int64> GetSourceFunctionIndexOfNode(const Node* node) const;

  // Returns the source function index of all nodes that map to 'map_target'.
  absl::StatusOr<std::set<int64>> GetSourceFunctionIndexesOfNodesMappedToNode(
      const Node* map_target) const;

  // Returns true if node_a and node_b are both map targets for nodes from a
  // common source function (or are themselves in the common function).
  absl::StatusOr<bool> NodeSourceFunctionsCollide(const Node* node_a,
                                                  const Node* node_b) const;

  // Returns a vector of Nodes to which the operands of the node
  // 'node' map. If node is owned by the integrated function, these are just
  // node's operands. If an operand does not yet have a mapping, the operand is
  // temporarily mapped to a new parameter(not yet implemented). Use of this
  // temporary will be replaced with the real mapping when it is set.
  absl::StatusOr<std::vector<Node*>> GetIntegratedOperands(
      const Node* node) const;

  // Returns true if 'node' is mapped to a node in the integrated function.
  bool HasMapping(const Node* node) const;

  // Returns true if all operands of 'node' are mapped to a node in the
  // integrated function.
  bool AllOperandsHaveMapping(const Node* node) const;

  // Returns true if other nodes map to 'node'
  bool IsMappingTarget(const Node* node) const;

  // Returns true if 'node' is in the integrated function.
  bool IntegrationFunctionOwnsNode(const Node* node) const {
    return function_.get() == node->function_base();
  }

  // Returns an estimate of the (gate count? area?) cost of a node.
  int64 GetNodeCost(const Node* node) const;

 private:
  IntegrationFunction(Package* package, const IntegrationOptions& options)
      : package_(package), integration_options_(options) {}

  // Helper function that implements the logic for merging nodes,
  // allowing for either the merge to be performed or for the cost
  // of the merge to be estimated.  A MergeNodesBackendResult struct is
  // returned. The 'can_merge' field indicates if the nodes can be merged. If
  // they can be merged, target_a and target_b point to resulting nodes that
  // represent the values of 'node_a' and 'node_b' in the integrated graph
  // (note: these will not necessarily point to the same merged node e.g. if the
  // merge node has a wider bitwidth than one of the original nodes, the target
  // pointer may instead point to a bitslice that takes in the wider node as an
  // operand). New muxes created or muxes modified by this call are placed in
  // 'changed_muxes'. Other nodes created by this call are placed in
  // 'other_added_nodes'.
  struct MergeNodesBackendResult {
    bool can_merge;
    // Separate targets for node_a and node_b because we may derive
    // map target nodes from the merged node;
    Node* target_a = nullptr;
    Node* target_b = nullptr;
    std::vector<UnifiedNode> changed_muxes;
    // We use a list rather than a vector here because
    // we will later want to remove elements in a (currently)
    // unknown order.  This would involve wastefule data copying
    // if we used a vector.
    std::list<Node*> other_added_nodes;
  };
  absl::StatusOr<MergeNodesBackendResult> MergeNodesBackend(const Node* node_a,
                                                            const Node* node_b);

  // For the integration function nodes node_a and node_b,
  // returns a single integration function node that combines the two
  // nodes. If a node combining node_a and node_b does not already
  // exist, a new mux and a per-mux select parameter are added.
  // If provided, the bool pointed to  by 'new_mux_added' will be
  // set to true if this call added a new mux.  Otherwise, false.
  absl::StatusOr<UnifiedNode> UnifyIntegrationNodesWithPerMuxSelect(
      Node* node_a, Node* node_b);

  // For the integration function nodes node_a and node_b,
  // returns a single integration function node that combines the two
  // nodes. If neither node_a or node_b is a mux added by a previous call,
  // a new mux is added whose select signal global_mux_select_param_.
  // If one of the nodes is such a mux, the other node is added as an
  // input to the mux.
  absl::StatusOr<UnifiedNode> UnifyIntegrationNodesWithGlobalMuxSelect(
      Node* node_a, Node* node_b);

  // Helper function for UnifyIntegrationNodesWithGlobalMuxSelect that handles
  // the case that neither input is a pre-existing mux.
  absl::StatusOr<UnifiedNode> UnifyIntegrationNodesWithGlobalMuxSelectArgIsMux(
      Node* node_a, Node* node_b);

  // Helper function for UnifyIntegrationNodesWithGlobalMuxSelect that handles
  // the cases that one of the input nodes is a pre-existing mux. The other
  // input is a node that will be added as a case(s) to the mux.
  absl::StatusOr<UnifiedNode> UnifyIntegrationNodesWithGlobalMuxSelectNoMuxArg(
      Node* mux, Node* case_node);

  // For a mux produced by UnifyIntegrationNodesWithPerMuxSelect, remove the mux
  // and the select parameter. Updates internal unification / mux book-keeping
  // accordingly. This function should only be called if the mux has no users
  // and mux's select signal is only used by the mux.
  absl::Status DeUnifyIntegrationNodesWithPerMuxSelect(Node* node);

  // For a mux produced by UnifyIntegrationNodesWithGlobalMuxSelect,
  // remove the most recently added cases. If no cases are in use after this,
  // then the mux is removed. Otherwise, the revert mux is returned. Also
  // updates internal unification / mux book-keeping.
  absl::StatusOr<Node*> DeUnifyIntegrationNodesWithGlobalMuxSelect(Node* node);

  // Replaces 'mux' with a new mux which is identical except that
  // the the cases at the indexes in 'source_index_to_case'
  // are replaced with the nodes specified.
  absl::StatusOr<Node*> ReplaceMuxCases(
      Node* mux_node,
      const absl::flat_hash_map<int64, Node*>& source_index_to_case);

  // Track mapping of original function nodes to integrated function nodes.
  absl::flat_hash_map<const Node*, Node*> original_node_to_integrated_node_map_;
  absl::flat_hash_map<const Node*, absl::flat_hash_set<const Node*>>
      integrated_node_to_original_nodes_map_;

  // Track which node-pairs have an associated mux.
  absl::flat_hash_map<std::pair<const Node*, const Node*>, Node*>
      node_pair_to_mux_;

  // If integration_options_ does not specify that each mux has a unique
  // select signal, this shared parameter is the select for all integration
  // muxes.
  Node* global_mux_select_param_ = nullptr;

  struct GlobalMuxMetadata {
    // Prefer to use sets rather than flat_hash_sets so that
    // the order of indexes reflects the order of the mux inputs.

    // Track which select arms map meaningful nodes for muxes
    // using the global_mux_select_param_ select signal.
    std::set<int64> occupied_case_indexes;

    // Track which select arms were most recently added for muxes
    // using the global_mux_select_param_ select signal. Note that
    // we only need to preserve enough history to call DeUnifyIntegrationNodes
    // after calling UnifyIntegrationNodes, with no other calls to
    // UnifyIntegrationNodes for a given mux inbetween. Further, there should
    // not be repeated calls to DeUnifyIntegrationNodes for a given node.
    std::set<int64> last_case_indexes_added;
  };
  // Track information about muxes that use global_mux_select_param_ as their
  // select signal.
  absl::flat_hash_map<Node*, GlobalMuxMetadata> global_mux_to_metadata_;

  // Source function in the integration package.
  std::vector<const Function*> source_functions_;

  // Maps each source function to a unique index.
  absl::flat_hash_map<const FunctionBase*, int64>
      source_function_base_to_index_;

  // Integrated function.
  std::unique_ptr<Function> function_;
  Package* package_;
  const IntegrationOptions integration_options_;
};

// A abstract class for an algorithm to integrate multiple functions.
template <class AlgorithmType>
class IntegrationAlgorithm {
 public:
  IntegrationAlgorithm(const IntegrationAlgorithm& other) = delete;
  void operator=(const IntegrationAlgorithm& other) = delete;

  // Returns a function that integrates the functions in source_functions.
  static absl::StatusOr<std::unique_ptr<IntegrationFunction>>
  IntegrateFunctions(Package* package,
                     absl::Span<const Function* const> source_functions,
                     const IntegrationOptions& options = IntegrationOptions());

 protected:
  IntegrationAlgorithm(Package* package,
                       absl::Span<const Function* const> source_functions,
                       const IntegrationOptions& options = IntegrationOptions())
      : package_(package), integration_options_(options) {
    source_functions_.reserve(source_functions.size());
    for (const auto* func : source_functions) {
      source_functions_.push_back(func);
    }
  }

  // Represents a possible modification to an integration function.
  enum class IntegrationMoveType { kInsert, kMerge };
  struct IntegrationMove {
    Node* node;
    IntegrationMoveType move_type;
    Node* merge_node = nullptr;
    int64 cost;
  };

  // Perform the modification to 'integration_function' described by 'move'.
  absl::StatusOr<std::vector<Node*>> ExecuteMove(
      IntegrationFunction* integration_function, const IntegrationMove& move);

  // Initialize any member fields.
  virtual absl::Status Initialize() = 0;

  // Returns a function that integrates the functions in source_functions_.
  // Runs after Initialize.
  virtual absl::StatusOr<std::unique_ptr<IntegrationFunction>> Run() = 0;

  // Create and return a new IntegrationFunction.
  absl::StatusOr<std::unique_ptr<IntegrationFunction>>
  NewIntegrationFunction() {
    return IntegrationFunction::MakeIntegrationFunctionWithParamTuples(
        package_, source_functions_, integration_options_);
  }

  // Get the IntegrationOptions::Algorithm value corresponding to
  // this algoirthm.
  virtual IntegrationOptions::Algorithm
  get_corresponding_algorithm_option() = 0;

  std::vector<const Function*> source_functions_;
  Package* package_;
  const IntegrationOptions integration_options_;
};

// A naive merging algorithm.  A node is eligible to be added to the
// integration function when all of its operands have already been added.
// At each step, adds the eligible node for which the cost of adding it to
// the function (either by inserting or merging with any integeration function
// node) is the lowest.
class BasicIntegrationAlgorithm
    : public IntegrationAlgorithm<BasicIntegrationAlgorithm> {
 public:
  BasicIntegrationAlgorithm(const BasicIntegrationAlgorithm& other) = delete;
  void operator=(const BasicIntegrationAlgorithm& other) = delete;

 private:
  // IntegrationAlgorithm::IntegrateFunctions needs to be able to call derived
  // class constructor.
  friend class IntegrationAlgorithm;

  BasicIntegrationAlgorithm(
      Package* package, absl::Span<const Function* const> source_functions,
      const IntegrationOptions& options = IntegrationOptions())
      : IntegrationAlgorithm(package, source_functions, options) {}

  // Represents a possible modification to the integration function.
  struct BasicIntegrationMove : IntegrationMove {
    std::list<Node*>::iterator node_itr;
  };

  // Make a BasicIntegrationMove for an insert.
  inline BasicIntegrationMove MakeInsertMove(
      std::list<Node*>::iterator node_itr, int64 cost) {
    return BasicIntegrationMove{{.node = *node_itr,
                                 .move_type = IntegrationMoveType::kInsert,
                                 .cost = cost},
                                .node_itr = node_itr};
  }

  // Make a BasicIntegrationMove for a merge.
  inline BasicIntegrationMove MakeMergeMove(std::list<Node*>::iterator node_itr,
                                            Node* merge_node, int64 cost) {
    return BasicIntegrationMove{{.node = *node_itr,
                                 .move_type = IntegrationMoveType::kMerge,
                                 .merge_node = merge_node,
                                 .cost = cost},
                                .node_itr = node_itr};
  }

  // Initialize member fields.
  absl::Status Initialize();

  // Returns a function that integrates the functions in source_functions_.
  // Runs after Initialize.
  absl::StatusOr<std::unique_ptr<IntegrationFunction>> Run();

  // Get the IntegrationOptions::Algorithm value corresponding to
  // this algoirthm.
  IntegrationOptions::Algorithm get_corresponding_algorithm_option() {
    return IntegrationOptions::Algorithm::kBasicIntegrationAlgorithm;
  }

  // Queue node for processing if all its operands are mapped
  // and node has not already been queued for processing.
  void EnqueueNodeIfReady(Node* node);

  // Track nodes for which all operands are already mapped and
  // are ready to be added to the integration_function_
  std::list<Node*> ready_nodes_;

  // Track all nodes that have ever been inserted into 'ready_nodes_'.
  absl::flat_hash_set<Node*> queued_nodes_;

  // Function combining the source functions.
  std::unique_ptr<IntegrationFunction> integration_function_;
};

// Class used to integrate separate functions into a combined, reprogrammable
// circuit that can be configured to have the same functionality as the
// input functions. The builder will attempt to construct the integrated
// function such that hardware common to the input functions is consolidated.
// Note that this is distinct from function inlining. With inlining, a function
// call is replaced by the body of the function that is called.  With function
// integration, we take separate functions that do not call each other and
// combine the hardware used to implement the functions.
class IntegrationBuilder {
 public:
  IntegrationBuilder(const IntegrationBuilder& other) = delete;
  void operator=(const IntegrationBuilder& other) = delete;

  // Creates an IntegrationBuilder and uses it to produce an integrated function
  // implementing all functions in source_functions_.
  static absl::StatusOr<std::unique_ptr<IntegrationBuilder>> Build(
      absl::Span<const Function* const> input_functions,
      const IntegrationOptions& options = IntegrationOptions());

  // Return functions to be integrated, in the integration package.
  absl::Span<Function*> source_functions() {
    return absl::Span<Function*>(source_functions_);
  }

  Package* package() { return package_.get(); }
  const IntegrationFunction* integrated_function() {
    return integrated_function_.get();
  }
  const IntegrationOptions* integration_options() {
    return &integration_options_;
  }

 private:
  IntegrationBuilder(absl::Span<const Function* const> input_functions,
                     const IntegrationOptions& options) {
    original_package_source_functions_.insert(
        original_package_source_functions_.end(), input_functions.begin(),
        input_functions.end());
    // TODO(jbaileyhandle): Make package name an optional argument.
    package_ = absl::make_unique<Package>("IntegrationPackage");
  }

  // Copy the source functions into a common package.
  absl::Status CopySourcesToIntegrationPackage();

  // Recursively copy a function into the common package_.
  absl::StatusOr<Function*> CloneFunctionRecursive(
      const Function* function,
      absl::flat_hash_map<const Function*, Function*>* call_remapping);

  // Set the integrated_function_.
  void set_integrated_function(
      std::unique_ptr<IntegrationFunction> integrated) {
    integrated_function_ = std::move(integrated);
  }

  // Uniquer to avoid function name collisions.
  NameUniquer function_name_uniquer_ = NameUniquer(/*separator=*/"__");

  // Options dictating how to integrate functions.
  const IntegrationOptions integration_options_;

  // Common package for to-be integrated functions
  // and integrated function.
  std::unique_ptr<Package> package_;

  // Function (and metadata) combining the source functions.
  std::unique_ptr<const IntegrationFunction> integrated_function_;

  // Functions to be integrated, in the integration package.
  std::vector<Function*> source_functions_;
  // Functions to be integrated, in their original packages.
  std::vector<const Function*> original_package_source_functions_;
};

}  // namespace xls

#endif  // XLS_INTEGRATOR_IR_INTEGRATOR_H_
