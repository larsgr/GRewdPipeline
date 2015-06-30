# Wrapper script
source( "${RSourceFile}" )

do.call( ${FunctionName}, args=as.list(commandArgs(trailingOnly = T)))
