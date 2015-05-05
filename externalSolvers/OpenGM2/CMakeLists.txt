 add_executable(dawid dawid_skene.cxx ${headers})
 if(WITH_HDF5)    
   target_link_libraries(dawid ${HDF5_LIBRARIES})  
 endif()
