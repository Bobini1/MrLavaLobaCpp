install(
    TARGETS MrLavaLoba_exe
    RUNTIME COMPONENT MrLavaLoba_Runtime
)

if(PROJECT_IS_TOP_LEVEL)
  include(CPack)
endif()
