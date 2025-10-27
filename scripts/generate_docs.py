import io

from contextlib import redirect_stdout
from SSDraw.cli import get_args


def generate_help(command, filename):
    f = io.StringIO()
    with redirect_stdout(f):
        try:
            get_args(argv=[command, "--help"])
        except SystemExit:
            pass

    help_output = f.getvalue()

    with open(f"docs/commands/{filename}.md", "w") as doc_file:
        doc_file.write(f"# SSDraw {command} Command\n\n")
        doc_file.write(
            f"This section documents the `ssdraw {command}` command.\n\n"
        )
        doc_file.write("```bash\n")
        doc_file.write(help_output)
        doc_file.write("```\n")


if __name__ == "__main__":
    generate_help("single", "single")
    generate_help("multi", "multi")

    print("CLI documentation updated in docs/commands/")
