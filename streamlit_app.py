import streamlit as st
import pandas as pd
import os
import io
import uuid
import base64

# Configuração da página
st.set_page_config(
    page_title="Consolidador de Genomas",
    page_icon="🧬",
    layout="wide"
)

# CSS personalizado
st.markdown("""
<style>
    .main {
        padding: 2rem;
    }
    .stAlert {
        border-radius: 10px;
    }
    .css-18e3th9 {
        padding-top: 1rem;
    }
    .css-1d391kg {
        padding: 1rem;
    }
    h1, h2, h3 {
        color: #2c3e50;
    }
    .stButton>button {
        background-color: #4a6fa5;
        color: white;
        border-radius: 8px;
        padding: 0.5rem 1rem;
        font-weight: bold;
    }
    .stButton>button:hover {
        background-color: #3a5a80;
    }
    .info-box {
        background-color: #e9f0f9;
        padding: 1rem;
        border-radius: 10px;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)

# Funções de processamento
def detectar_separador(conteudo):
    """Detecta automaticamente o separador do arquivo CSV"""
    # Obtém a primeira linha do conteúdo
    primeira_linha = conteudo.splitlines()[0] if conteudo.splitlines() else ""
    
    # Conta a quantidade de cada possível separador
    virgula = primeira_linha.count(',')
    ponto_virgula = primeira_linha.count(';')
    
    # Retorna o separador mais frequente
    if ponto_virgula > virgula:
        return ';'
    return ','

def consolidar_genomas_csv(df):
    """Consolida os dados de genomas no DataFrame"""
    # Verifica se 'record_id' existe no arquivo
    if 'record_id' not in df.columns:
        st.error("Erro: Coluna 'record_id' não encontrada no arquivo!")
        return None, 0, 0
    
    # Obter a lista de colunas do arquivo
    colunas = list(df.columns)
    
    # Registros originais
    registros_originais = len(df)
    
    # Obter a lista de record_ids únicos
    record_ids_unicos = df['record_id'].unique()
    
    # Barra de progresso
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    # Para cada record_id, consolidar as informações
    resultados = []
    total = len(record_ids_unicos)
    
    for i, record_id in enumerate(record_ids_unicos):
        # Atualizar barra de progresso
        progress = int((i + 1) / total * 100)
        progress_bar.progress(progress)
        status_text.text(f"Processando... {progress}% concluído ({i+1}/{total})")
        
        # Filtrar as linhas para o record_id atual
        linhas_amostra = df[df['record_id'] == record_id]
        
        # Criar um dicionário para armazenar os valores consolidados
        amostra_consolidada = {'record_id': record_id}
        
        # Inicializar os outros campos como vazios
        for coluna in colunas:
            if coluna != 'record_id':
                amostra_consolidada[coluna] = ''
        
        # Para cada coluna, preencher com valores não vazios
        for coluna in colunas:
            if coluna != 'record_id':
                # Obter valores não nulos/vazios
                valores = linhas_amostra[coluna].dropna().astype(str)
                valores = valores[valores != '']
                
                # Se houver valores, use o primeiro
                if len(valores) > 0:
                    amostra_consolidada[coluna] = valores.iloc[0]
        
        resultados.append(amostra_consolidada)
    
    # Limpar a barra de progresso e o texto de status
    progress_bar.empty()
    status_text.empty()
    
    # Converter a lista de dicionários em um DataFrame
    df_consolidado = pd.DataFrame(resultados)
    
    # Garantir que a ordem das colunas seja mantida
    df_consolidado = df_consolidado[colunas]
    
    return df_consolidado, registros_originais, len(df_consolidado)

def get_csv_download_link(df, filename="dados_consolidados.csv", text="Baixar CSV"):
    """Gera um link para download do DataFrame como CSV"""
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}" class="download-button">{text}</a>'
    return href

# Interface do aplicativo
st.title("Consolidador de Genomas 🧬")
st.markdown("### Ferramenta para consolidar registros de genomas com base no 'record_id'")

# Caixa de informações
with st.expander("ℹ️ Instruções e Informações", expanded=True):
    st.markdown("""
    Esta ferramenta consolida múltiplas linhas em um arquivo CSV que compartilham o mesmo 'record_id'.
    
    **Como funciona:**
    1. Faça upload de um arquivo CSV contendo dados genômicos
    2. O arquivo deve conter uma coluna chamada 'record_id'
    3. Para cada 'record_id', será mantido apenas um registro, preservando a primeira ocorrência de cada valor
    4. O resultado será exibido na tela e disponibilizado para download
    
    **Dicas:**
    - O separador (vírgula ou ponto-e-vírgula) é detectado automaticamente
    - Arquivos muito grandes podem levar mais tempo para processar
    - Todos os dados são processados apenas no seu navegador, nada é armazenado em servidores externos
    """)

# Upload de arquivo
uploaded_file = st.file_uploader("Selecione um arquivo CSV", type=["csv"])

if uploaded_file is not None:
    try:
        # Ler o conteúdo do arquivo
        conteudo = uploaded_file.getvalue().decode('utf-8')
        
        # Detectar o separador
        separador = detectar_separador(conteudo)
        
        # Ler o arquivo CSV
        df = pd.read_csv(io.StringIO(conteudo), sep=separador)
        
        # Mostrar informações básicas
        st.write(f"**Arquivo carregado:** {uploaded_file.name}")
        st.write(f"**Separador detectado:** '{separador}'")
        st.write(f"**Número de linhas:** {len(df)}")
        st.write(f"**Número de colunas:** {len(df.columns)}")
        
        # Botão de processamento
        if st.button("🔄 Processar Arquivo"):
            with st.spinner("Processando o arquivo... Por favor, aguarde."):
                # Processar o arquivo
                df_consolidado, registros_originais, registros_consolidados = consolidar_genomas_csv(df)
                
                if df_consolidado is not None:
                    # Calcular porcentagem de redução
                    reducao_percentual = round(((registros_originais - registros_consolidados) / registros_originais * 100), 2) if registros_originais > 0 else 0
                    
                    # Mostrar resultados em colunas
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Registros Originais", registros_originais)
                    with col2:
                        st.metric("Registros Consolidados", registros_consolidados)
                    with col3:
                        st.metric("Redução", f"{reducao_percentual}%")
                    
                    # Mostrar preview
                    st.subheader("Preview dos Dados Consolidados")
                    st.dataframe(df_consolidado.head(10))
                    
                    # Gerar link de download
                    st.markdown("### Download")
                    nome_arquivo = f"{os.path.splitext(uploaded_file.name)[0]}_consolidado.csv"
                    st.markdown(get_csv_download_link(df_consolidado, nome_arquivo, "📥 Baixar Arquivo Consolidado"), unsafe_allow_html=True)
                    
                    # Mostrar informações adicionais
                    st.info(f"""
                    **Informações adicionais:**
                    - O arquivo consolidado mantém todas as {len(df.columns)} colunas do original
                    - Apenas os primeiros valores não vazios de cada registro foram preservados
                    - O arquivo gerado ({nome_arquivo}) pode ser aberto em qualquer programa que leia CSV
                    """)
    
    except Exception as e:
        st.error(f"Ocorreu um erro ao processar o arquivo: {str(e)}")
        st.info("""
        **Possíveis soluções:**
        - Verifique se o arquivo CSV contém uma coluna chamada 'record_id'
        - Certifique-se de que o arquivo está no formato CSV correto
        - Verifique se o arquivo não está corrompido
        - Tente novamente com outro arquivo CSV
        """)

# Rodapé
st.markdown("---")
st.markdown("Versão 1.0 - Desenvolvido para processamento de dados genômicos")
