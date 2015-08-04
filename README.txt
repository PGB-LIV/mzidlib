=========================================================================================================================
=                                                                     													=
=                                                                     													=
=                MzidLib - A library of routines for manipulating data in the PSI's mzIdentML standard.                 =
=                                                                     													=
=                                                                     													=
=========================================================================================================================



Description

  MzidLib - A library of routines for manipulating data in the PSI's mzIdentML standard. 

How to use in maven
  
  mzidlib artifact is currently using Googlecode host as its maven
  repository. To use it, add the following to your POM:
  
  First, add MzidLib dependency (check the latest version):
  
  <dependency>
     <groupId>MzidLib</groupId>
     <artifactId>mzidlib</artifactId>
     <version>1.6</version>
  </dependency>

  Second, add wagon-svn extension (current version 1.12) in <build> tag:
  
  <extensions>
     <extension>
        <groupId>org.jvnet.wagon-svn</groupId>
        <artifactId>wagon-svn</artifactId>
        <version>1.12</version>
     </extension>
  </extensions>

  Finally, add plugin repository and mzidlib repository:
  
  <pluginRepositories>
     <pluginRepository>
        <id>maven2-repository.dev.java.net</id>
        <name>Java.net Repository for Maven</name>
        <url>http://download.java.net/maven/2/</url>
     </pluginRepository>
  </pluginRepositories>
  <repositories>
     <repository>
        <id>mzidentml-lib-maven-repo</id>
        <name>Maven Repository for mzidlib release</name>
        <url>http://mzidentml-lib.googlecode.com/svn/maven/repo</url>
     </repository>
  </repositories>


=============================================================================================================
Type java -jar mzidentml-lib.jar to get a description of the tools and command line parameters for running the mzidLibrary. 

Please download example files from https://mzidentml-lib.googlecode.com/files/example_files.zip and place these in a directory alongside test_outputs (see run_examples.bat for location).
   
  