// ---------------------------------------------------------------------
//
// Copyright (C) 2026 by Luca Heltai
//
// This file is part of the deal.II Testbench for NMOPT Laboratories,
// deal.II library.
//
// The deal.II Testbench for NMOPT Laboratories is free software; you can use
// it, redistribute it, and/or modify it under the terms of the Apache-2.0
// License WITH LLVM-exception as published by the Free Software Foundation;
// either version 3.0 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at the top
// level of the deal.II Testbench for NMOPT Laboratories distribution.
//
// ---------------------------------------------------------------------

#include "inverse_poisson_kkt.h"

#include <deal.II/base/parameter_acceptor.h>

#include <iostream>

using namespace dealii;


int
main(int argc, char **argv)
{
  try
    {
      InversePoissonKKT<DEAL_DIMENSION> inverse_problem;

      const std::string parameter_file =
        (argc > 1 ? argv[1] : "parameters/inverse_poisson_kkt.prm");

      ParameterAcceptor::initialize(parameter_file, "last_used_parameters.prm");
      inverse_problem.run();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
