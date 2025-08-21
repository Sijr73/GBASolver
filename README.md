
# GBA Solver

## Overview

**GBA Solver**  is an intuitive web application for Growth Balance Analysis (GBA) of self-replicating cell models. Built with R/Shiny, it features a user-friendly, spreadsheet-style interface that enables researchers to model cellular resource allocation under nonlinear kinetic rate laws. GBA Solver streamlines model construction by integrating enzyme kinetics and parameter data from the BRENDA database, and generating dynamic visualizations of metabolic pathways. By lowering computational barriers, GBA Solver makes nonlinear cellular modeling more accessible to the broader scientific community.

🔗 **Try GBA Solver online:** [GBA Solver Web Application](https://cellgrowthsim.com/) 

This repository contains all the necessary source code and Docker configuration files for building and running GBA Solver locally or on a server.

1. **Source Code** for the GBA Solver Shiny application (in the `GBApp` directory)  
2. **Dockerfile** for building and running the GBA Solver container locally or on a server

---

##  Running GBA Solver locally via Docker

These instructions let you build and run the GBA Solver Docker container on your local machine.

1. **Install Docker**  
   - [Docker Desktop](https://www.docker.com/) for macOS, Windows, or Linux.

2. **Clone or Download this Repository**

   ```bash
   git clone https://github.com/Sijr73/CellGrowthSimulator.git
   cd GBApp

3. **Build the Docker Image**

   
   ```bash
   docker build -t gbapp .
   ```

4. **Run the Docker Container**

   ```bash
   
   docker run -d --rm -p 3838:3838 --name gbapp gbapp
   ```
   `-d`: run the container in detached mode

   `-p 3838:3838`: map port 3838 in the container to port 3838 on your machine

5. **Open GBA Solver**

   - Go to http://localhost:3838 in your web browser.
   - You should see the GBA Solver home page.
## Contact
If you have any questions or run into issues, please open an issue in the repository or contact the repository maintainer.
