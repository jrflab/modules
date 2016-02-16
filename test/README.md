# Test
## Saba test
```
./test.sh
```
## Local test
- Mount saba over sshfs
- create symbolic link from local ~/share to sshfs mounted share dir
- create anaconda environment
    ```
    conda create -p test-env --file ../conda_env/jrflab_modules_env.txt
    ```
- run test
    ```
    export USE_CLUSTER=false JRFLAB_MODULES_ENV=fullpathtotest-env && ./test.sh 
    ```
