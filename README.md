
# Cell Growth Simulator

## Overview

**Cell Growth Simulator**  is an intuitive web application for Growth Balance Analysis (GBA) of self-replicating cell models. Built with R/Shiny, it features a user-friendly, spreadsheet-style interface that enables researchers to model cellular resource allocation under nonlinear kinetic rate laws. Cell Growth Simulator streamlines model construction by integrating enzyme kinetics and parameter data from the BRENDA database, and generating dynamic visualizations of metabolic pathways. By lowering computational barriers, Cell Growth Simulator makes nonlinear cellular modeling more accessible to the broader scientific community.

🔗 **Try Cell Growth Simulator online:** [Cell Growth Simulator Web Application](https://cellgrowthsim.com/) 

This repository contains all the necessary source code and Docker configuration files for building and running Cell Growth Simulator locally or on a server.

1. **Source Code** for the Cell Growth Simulator Shiny application (in the `GBApp` directory)  
2. **Dockerfile** for building and running the Cell Growth Simulator container locally or on a server

---

##  Running Cell Growth Simulator locally via Docker

These instructions let you build and run the Cell Growth Simulator Docker container on your local machine.

1. **Install Docker**  
   - [Docker Desktop](https://www.docker.com/) for macOS, Windows, or Linux.

2. **Clone or Download this Repository**

   ```bash
   git clone https://github.com/Sijr73/GBApp.git
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

5. **Open GBApp**

   - Go to http://localhost:3838 in your web browser.
   - You should see the GBApp home page.
## Contact
If you have any questions or run into issues, please open an issue in the repository or contact the repository maintainer.
