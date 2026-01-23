# Model Predictive Control

- [Model Predictive Control](#model-predictive-control)
  - [Overview](#overview)
  - [Setup](#setup)
    - [Dependencies](#dependencies)
  - [Docker](#docker)
  - [Reference](#reference)


## Overview

This project implements the different variants of MPC:
- Linear MPC
- Linear Time-Varying MPC


## Setup

### Dependencies
- Eigen
    ```sh
    sudo apt install libeigen3-dev
    ```
- SDL2
    ```sh
    sudo apt install libsdl2-dev libsdl2-ttf-dev
    ```

## Docker
- Add: `xhost +local:docker`
- Build docker image
    ```sh
    docker compose build
    ```
- Run and enter the docker container
    ```sh
    docker compose run --rm model_predictive_control bash
    ```


## Reference

