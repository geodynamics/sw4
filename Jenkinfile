#!groovy

pipeline {
  agent {
    docker {
      image 'geodynamics/sw4-buildenv-bionic:latest'
      alwaysPull true
    }
  }

  options {
    timeout(time: 4, unit: 'HOURS')
  }

  stages {
    stage('Build') {
      steps {
        sh 'mkdir build'
        sh '''
          cd build

          # MPICH has issues when there are more processes than hardware threads
          # See https://wiki.mpich.org/mpich/index.php/Frequently_Asked_Questions#Q:_Why_does_my_MPI_program_run_much_slower_when_I_use_more_processes.3F
          cmake \
            -D MPI_NUM_TEST_PROCS='4' \
            -D TESTING_LEVEL='1' \
            ..
        '''
        sh '''
          cd build
          make
        '''
      }
    }

    stage('Test') {
      steps {
        sh '''
          cd build
          ctest --no-compress-output -T Test
        '''
      }
      post {
        always {
          xunit testTimeMargin: '3000',
          thresholdMode: 1,
          thresholds: [failed(), skipped()],
          tools: [CTest(pattern: 'build/Testing/**/*.xml')]
        }
      }
    }
  }

  post { always { cleanWs() } }
}
