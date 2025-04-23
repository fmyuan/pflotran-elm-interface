if [ -n "$ARTIFACT_DIR" ] && 
   [ $(grep -c "failed" "$ARTIFACT_DIR/status") -ne 0 ]; then
  echo 'A prior build/test failed. Skipping current.'
  exit 0
fi
