#!/bin/bash

EXEC=./exato4
INPUT_DIR=./output_folder

# ---------- Checks ----------
if [ ! -x "$EXEC" ]; then
  echo "Erro: executável $EXEC não encontrado ou sem permissão."
  exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
  echo "Erro: diretório $INPUT_DIR não existe."
  exit 1
fi

echo "Iniciando execução em ordem alfabética..."
echo

# ---------- Loop ordenado ----------
find "$INPUT_DIR" -maxdepth 1 -type f \
  | sort -V \
  | while IFS= read -r file; do

      base=$(basename "$file")

      # ignora arquivos ocultos
      [[ "$base" == .* ]] && continue

      echo "======================================"
      echo "Rodando instância: $base"
      echo "Comando: $EXEC $file"
      echo "======================================"

      "$EXEC" "$file"
      status=$?

      if [ $status -ne 0 ]; then
        echo "⚠️  Erro ao processar $base (exit code $status)"
      fi

      echo
    done

echo "✅ Todas as instâncias foram processadas."
