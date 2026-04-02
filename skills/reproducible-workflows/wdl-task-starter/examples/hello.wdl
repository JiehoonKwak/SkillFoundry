version 1.0

task write_greeting {
  input {
    String name
  }

  command <<<
    printf "Hello, ~{name}\n" > greeting.txt
  >>>

  output {
    File greeting_file = "greeting.txt"
    String greeting_text = read_string("greeting.txt")
  }
}

workflow hello_workflow {
  input {
    String name
  }

  call write_greeting {
    input:
      name = name
  }

  output {
    File greeting_file = write_greeting.greeting_file
    String greeting_text = write_greeting.greeting_text
  }
}
