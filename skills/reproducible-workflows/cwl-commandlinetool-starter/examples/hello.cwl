cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - python3
  - -c
  - |
    from pathlib import Path
    import sys
    Path("greeting.txt").write_text(sys.argv[1] + "\n", encoding="utf-8")
inputs:
  message:
    type: string
    inputBinding:
      position: 1
outputs:
  greeting_file:
    type: File
    outputBinding:
      glob: greeting.txt
