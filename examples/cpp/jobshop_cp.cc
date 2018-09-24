// Copyright 2018 Google LLC
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "ortools/base/logging.h"
#include "ortools/constraint_solver/constraint_solver.h"

// Solve a job shop problem:

namespace operations_research {

  void SolveJobShopExample() {
    // Instantiate the solver.
    Solver solver("JobShopExample");
    std::array<int, 3> machines;
    std::iota(std::begin(machines), std::end(machines), 0);

    std::array<int, 3> jobs;
    std::iota(std::begin(jobs), std::end(jobs), 0);

    // Array of Jobs
    using MachineIndex = int;
    using ProcessingTime = int;
    using Task = std::pair<MachineIndex, ProcessingTime>;
    using Job = std::vector<Task>;
    std::array<Job, jobs.size()> jobs_list = {
      {{0, 3}, {1, 2}, {2, 2}},
      {{0, 2}, {2, 1}, {1, 4}},
      {{1, 4}, {2, 3}}
    };

    // Computes horizon.
    ProcessingTime horizon = 0;
    for (const Job& job: jobs_list) {
      for (const Task& task: job) {
        horizon += task.second;
      }
    }

    // Creates tasks.
    std::array<std::vector<IntVar*>, jobs.size()> tasks_matrix;
    for (const auto job: jobs) {
      for (std::size_t task=0; task < jobs_list[job].size(); ++task) {
				std::ostringstream oss;
				oss << "Tasks(job: " << job << ", task: " << task << ")";
        tasks_matrix[job].push_back(
            solver.MakeFixedDurationIntervalVar(
              0,
              horizon,
              jobs[job][task].second,
              False,
              oss.str()));
      }
    }

    // Creates sequence variables and add disjunctive constraints.
    std::vector<SequenceVar*> all_sequences;
    std::vector<IntVar*> all_machines_jobs;
    for (const auto machine: machines) {
      std::vector<IntVar*> machines_jobs;
      for (const auto& j: jobs) {
        const Job& job = jobs_list[j];
        for (int task=0; task < job.size(); ++task) {
          if (job[task].fisrt == machine) machines_jobs.push_back(tasks_matrix[j, k]);
        }
      disj = solver.DisjunctiveConstraint(machines_jobs, 'machine %i' % i)
                  all_sequences.append(disj.SequenceVar())
                  solver.Add(disj)

                  // Add conjunctive contraints.
                  for i in all_jobs:
                    for j in range(0, len(machines[i]) - 1):
                      solver.Add(all_tasks[(i, j + 1)].StartsAfterEnd(all_tasks[(i, j)]))

                      // Set the objective.
                      obj_var = solver.Max([all_tasks[(i, len(machines[i])-1)].EndExpr()
                          for i in all_jobs])
                        objective_monitor = solver.Minimize(obj_var, 1)
                        // Create search phases.
                        sequence_phase = solver.Phase([all_sequences[i] for i in all_machines],
                            solver.SEQUENCE_DEFAULT)
                        vars_phase = solver.Phase([obj_var],
                            solver.CHOOSE_FIRST_UNBOUND,
                            solver.ASSIGN_MIN_VALUE)
                        main_phase = solver.Compose([sequence_phase, vars_phase])
                        // Create the solution collector.
                        collector = solver.LastSolutionCollector()

                        // Add the interesting variables to the SolutionCollector.
                        collector.Add(all_sequences)
                        collector.AddObjective(obj_var)

                        for i in all_machines:
                          sequence = all_sequences[i];
    sequence_count = sequence.Size();
    for j in range(0, sequence_count):
      t = sequence.Interval(j)
      collector.Add(t.StartExpr().Var())
      collector.Add(t.EndExpr().Var())
      // Solve the problem.
      disp_col_width = 10
      if solver.Solve(main_phase, [objective_monitor, collector]):
        print("\nOptimal Schedule Length:", collector.ObjectiveValue(0), "\n")
        sol_line = ""
        sol_line_tasks = ""
        print("Optimal Schedule", "\n")

        for i in all_machines:
          seq = all_sequences[i]
            sol_line += "Machine " + str(i) + ": "
            sol_line_tasks += "Machine " + str(i) + ": "
            sequence = collector.ForwardSequence(0, seq)
            seq_size = len(sequence)

            for j in range(0, seq_size):
              t = seq.Interval(sequence[j]);
    // Add spaces to output to align columns.
    sol_line_tasks +=  t.Name() + " " * (disp_col_width - len(t.Name()))

      for j in range(0, seq_size):
        t = seq.Interval(sequence[j]);
    sol_tmp = "[" + str(collector.Value(0, t.StartExpr().Var())) + ","
      sol_tmp += str(collector.Value(0, t.EndExpr().Var())) + "] "
      // Add spaces to output to align columns.
      sol_line += sol_tmp + " " * (disp_col_width - len(sol_tmp))

      sol_line += "\n"
      sol_line_tasks += "\n"

      print(sol_line_tasks)
      print("Time Intervals for Tasks\n")
      print(sol_line)
  }
}  // namespace operations_research

int main(int argc, char **argv) {
  google::InitGoogleLogging(argv[0]);
  FLAGS_logtostderr = 1;
  operations_research::SolveJobShopExample();
  return 0;
}

