import streamlit as st
import pandas as pd
import os
import io
import uuid
import base64
import plotly.express as px
import plotly.graph_objects as go
from collections import defaultdict

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
    .variant-card {
        border: 1px solid #ddd;
        border-radius: 8px;
        padding: 1rem;
        margin: 0.5rem 0;
        background-color: #f8f9fa;
    }
    .pathogenic {
        background-color: #ffebee;
        border-left: 4px solid #f44336;
    }
    .likely-pathogenic {
        background-color: #fff3e0;
        border-left: 4px solid #ff9800;
    }
    .benign {
        background-color: #e8f5e9;
        border-left: 4px solid #4caf50;
    }
    .likely-benign {
        background-color: #e0f2f1;
        border-left: 4px solid #009688;
    }
    .vus {
        background-color: #fff9c4;
        border-left: 4px solid #ffc107;
    }
    .conflicting {
        background-color: #f3e5f5;
        border-left: 4px solid #9c27b0;
    }
    .gene-header {
        font-weight: bold;
        color: #1976d2;
        font-size: 1.1em;
        margin-top: 1rem;
    }
    .case-header {
        font-weight: bold;
        color: #424242;
        font-size: 1.3em;
        margin-top: 1.5rem;
        padding: 0.5rem;
        background-color: #e3f2fd;
        border-radius: 5px;
    }
    .significance-header {
        font-weight: bold;
        color: #ffffff;
        font-size: 1.5em;
        margin-top: 2rem;
        padding: 1rem;
        border-radius: 8px;
        text-align: center;
    }
</style>
""", unsafe_allow_html=True)

# Funções de processamento
def detectar_separador(conteudo):
    """Detecta automaticamente o separador do arquivo CSV"""
    primeira_linha = conteudo.splitlines()[0] if conteudo.splitlines() else ""

    virgula = primeira_linha.count(',')
    ponto_virgula = primeira_linha.count(';')

    if ponto_virgula > virgula:
        return ';'
    return ','

def consolidar_genomas_csv(df):
    """Consolida os dados de genomas no DataFrame"""
    if 'record_id' not in df.columns:
        st.error("Erro: Coluna 'record_id' não encontrada no arquivo!")
        return None, 0, 0

    colunas = list(df.columns)
    registros_originais = len(df)
    record_ids_unicos = df['record_id'].unique()

    progress_bar = st.progress(0)
    status_text = st.empty()

    resultados = []
    total = len(record_ids_unicos)

    for i, record_id in enumerate(record_ids_unicos):
        progress = int((i + 1) / total * 100)
        progress_bar.progress(progress)
        status_text.text(f"Processando... {progress}% concluído ({i+1}/{total})")

        linhas_amostra = df[df['record_id'] == record_id]
        amostra_consolidada = {'record_id': record_id}

        for coluna in colunas:
            if coluna != 'record_id':
                amostra_consolidada[coluna] = ''

        for coluna in colunas:
            if coluna != 'record_id':
                valores = linhas_amostra[coluna].dropna().astype(str)
                valores = valores[valores != '']

                if len(valores) > 0:
                    amostra_consolidada[coluna] = valores.iloc[0]

        resultados.append(amostra_consolidada)

    progress_bar.empty()
    status_text.empty()

    df_consolidado = pd.DataFrame(resultados)
    df_consolidado = df_consolidado[colunas]

    return df_consolidado, registros_originais, len(df_consolidado)

def get_clnsig_color(clnsig):
    """Retorna a cor e classe CSS para cada tipo de significância clínica"""
    clnsig_lower = str(clnsig).lower()

    if 'pathogenic' in clnsig_lower and 'likely' not in clnsig_lower:
        return '#f44336', 'pathogenic', 'Pathogenic'
    elif 'likely_pathogenic' in clnsig_lower or 'likely pathogenic' in clnsig_lower:
        return '#ff9800', 'likely-pathogenic', 'Likely Pathogenic'
    elif 'benign' in clnsig_lower and 'likely' not in clnsig_lower:
        return '#4caf50', 'benign', 'Benign'
    elif 'likely_benign' in clnsig_lower or 'likely benign' in clnsig_lower:
        return '#009688', 'likely-benign', 'Likely Benign'
    elif 'uncertain' in clnsig_lower or 'vus' in clnsig_lower:
        return '#ffc107', 'vus', 'Uncertain Significance (VUS)'
    elif 'conflicting' in clnsig_lower:
        return '#9c27b0', 'conflicting', 'Conflicting Interpretations'
    else:
        return '#9e9e9e', 'vus', clnsig

def parse_genomic_data(df):
    """Processa dados genômicos do Excel"""
    required_columns = ['Case_Record_ID', 'Gene', 'CHR', 'POS', 'rsID',
                       'REF', 'ALT', 'Genotype', 'Zygosity', 'CLNSIG_Resolved']

    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        st.error(f"Colunas obrigatórias ausentes: {', '.join(missing_columns)}")
        return None

    # Limpar e normalizar dados
    df = df.copy()
    df['CLNSIG_Resolved'] = df['CLNSIG_Resolved'].fillna('Unknown')
    df['CLNDN'] = df['CLNDN'].fillna('Not provided') if 'CLNDN' in df.columns else 'Not provided'

    return df

def generate_genomic_report(df):
    """Gera relatório compreensivo categorizado por CLNSIG_Resolved e depois por casos"""

    # Agrupar por CLNSIG_Resolved
    clnsig_groups = df.groupby('CLNSIG_Resolved')

    # Estatísticas gerais
    st.markdown("## 📊 Resumo Estatístico")

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Total de Variantes", len(df))
    with col2:
        st.metric("Casos Únicos", df['Case_Record_ID'].nunique())
    with col3:
        st.metric("Genes Únicos", df['Gene'].nunique())
    with col4:
        st.metric("Categorias de Significância", df['CLNSIG_Resolved'].nunique())

    # Gráfico de distribuição de variantes por significância clínica
    st.markdown("### Distribuição por Significância Clínica")
    clnsig_counts = df['CLNSIG_Resolved'].value_counts()

    colors_map = {}
    for clnsig in clnsig_counts.index:
        color, _, _ = get_clnsig_color(clnsig)
        colors_map[clnsig] = color

    fig = go.Figure(data=[go.Bar(
        x=clnsig_counts.index,
        y=clnsig_counts.values,
        marker_color=[colors_map[x] for x in clnsig_counts.index],
        text=clnsig_counts.values,
        textposition='auto',
    )])
    fig.update_layout(
        title="Variantes por Significância Clínica",
        xaxis_title="Significância Clínica",
        yaxis_title="Número de Variantes",
        height=400
    )
    st.plotly_chart(fig, use_container_width=True)

    # Gráfico de genes mais frequentes
    st.markdown("### Top 20 Genes com Mais Variantes")
    gene_counts = df['Gene'].value_counts().head(20)
    fig2 = go.Figure(data=[go.Bar(
        x=gene_counts.index,
        y=gene_counts.values,
        marker_color='#4a6fa5',
        text=gene_counts.values,
        textposition='auto',
    )])
    fig2.update_layout(
        title="Genes Mais Frequentes",
        xaxis_title="Gene",
        yaxis_title="Número de Variantes",
        height=400
    )
    st.plotly_chart(fig2, use_container_width=True)

    # Relatório detalhado por CLNSIG_Resolved
    st.markdown("---")
    st.markdown("## 🧬 Relatório Detalhado por Significância Clínica")

    # Ordenar grupos por prioridade clínica
    priority_order = ['Pathogenic', 'Likely_pathogenic', 'Uncertain_significance',
                     'Likely_benign', 'Benign', 'Conflicting']

    for clnsig_name, group in sorted(clnsig_groups,
                                     key=lambda x: priority_order.index(x[0]) if x[0] in priority_order else 999):
        color, css_class, display_name = get_clnsig_color(clnsig_name)

        # Header de significância
        st.markdown(f"""
        <div class="significance-header" style="background-color: {color};">
            {display_name} ({len(group)} variantes)
        </div>
        """, unsafe_allow_html=True)

        # Agrupar por casos dentro de cada categoria de significância
        case_groups = group.groupby('Case_Record_ID')

        for case_id, case_data in case_groups:
            # Header do caso
            st.markdown(f"""
            <div class="case-header">
                📋 Caso: {case_id} ({len(case_data)} variantes)
            </div>
            """, unsafe_allow_html=True)

            # Agrupar por gene dentro de cada caso
            gene_groups = case_data.groupby('Gene')

            for gene, gene_data in gene_groups:
                st.markdown(f'<div class="gene-header">🧬 Gene: {gene}</div>',
                          unsafe_allow_html=True)

                # Mostrar cada variante
                for idx, variant in gene_data.iterrows():
                    zygosity_emoji = "🔴🔴" if variant['Zygosity'] == 'Homozygous' else "🔴⚪"

                    variant_html = f"""
                    <div class="variant-card {css_class}">
                        <strong>Posição:</strong> Chr{variant['CHR']}:{variant['POS']} |
                        <strong>rsID:</strong> {variant.get('rsID', 'N/A')} |
                        <strong>Zigosidade:</strong> {variant['Zygosity']} {zygosity_emoji}
                        <br>
                        <strong>Referência:</strong> {variant['REF']} → <strong>Alternativa:</strong> {variant['ALT']} |
                        <strong>Genótipo:</strong> {variant['Genotype']}
                        <br>
                        <strong>Condição Clínica:</strong> {variant.get('CLNDN', 'Not provided')}
                    """

                    if 'CLNSIGCONF' in variant and pd.notna(variant['CLNSIGCONF']) and variant['CLNSIGCONF'] != '':
                        variant_html += f"<br><strong>Conflito de Interpretação:</strong> {variant['CLNSIGCONF']}"

                    variant_html += "</div>"

                    st.markdown(variant_html, unsafe_allow_html=True)

    # Tabela de exportação
    st.markdown("---")
    st.markdown("## 📥 Exportar Dados")

    # Preparar dados para exportação
    export_df = df.copy()
    export_df = export_df.sort_values(['CLNSIG_Resolved', 'Case_Record_ID', 'Gene'])

    # Botão de download
    csv = export_df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="genomic_report.csv" style="text-decoration: none;"><button style="background-color: #4a6fa5; color: white; padding: 0.5rem 1rem; border: none; border-radius: 8px; cursor: pointer; font-weight: bold;">📥 Baixar Relatório Completo (CSV)</button></a>'
    st.markdown(href, unsafe_allow_html=True)

    # Excel download
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        export_df.to_excel(writer, index=False, sheet_name='Genomic Report')
    excel_data = output.getvalue()
    b64_excel = base64.b64encode(excel_data).decode()
    href_excel = f'<a href="data:application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;base64,{b64_excel}" download="genomic_report.xlsx" style="text-decoration: none;"><button style="background-color: #2e7d32; color: white; padding: 0.5rem 1rem; border: none; border-radius: 8px; cursor: pointer; font-weight: bold; margin-left: 1rem;">📥 Baixar Relatório Completo (Excel)</button></a>'
    st.markdown(href_excel, unsafe_allow_html=True)

    # Mostrar preview da tabela
    st.markdown("### Preview da Tabela Completa")
    st.dataframe(export_df, use_container_width=True)

def get_csv_download_link(df, filename="dados_consolidados.csv", text="Baixar CSV"):
    """Gera um link para download do DataFrame como CSV"""
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="{filename}" class="download-button">{text}</a>'
    return href

# Interface do aplicativo
st.title("Consolidador de Genomas 🧬")
st.markdown("### Ferramenta para consolidar registros de genomas e gerar relatórios de variantes genômicas")

# Seletor de modo
mode = st.radio(
    "Selecione o modo de operação:",
    ["Consolidar CSV (por record_id)", "Relatório de Variantes Genômicas (Excel)"],
    horizontal=True
)

if mode == "Consolidar CSV (por record_id)":
    # Modo original - Consolidação de CSV
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
        - Todos os dados são processados apenas no seu navegador
        """)

    uploaded_file = st.file_uploader("Selecione um arquivo CSV", type=["csv"])

    if uploaded_file is not None:
        try:
            conteudo = uploaded_file.getvalue().decode('utf-8')
            separador = detectar_separador(conteudo)
            df = pd.read_csv(io.StringIO(conteudo), sep=separador)

            st.write(f"**Arquivo carregado:** {uploaded_file.name}")
            st.write(f"**Separador detectado:** '{separador}'")
            st.write(f"**Número de linhas:** {len(df)}")
            st.write(f"**Número de colunas:** {len(df.columns)}")

            if st.button("🔄 Processar Arquivo"):
                with st.spinner("Processando o arquivo... Por favor, aguarde."):
                    df_consolidado, registros_originais, registros_consolidados = consolidar_genomas_csv(df)

                    if df_consolidado is not None:
                        reducao_percentual = round(((registros_originais - registros_consolidados) / registros_originais * 100), 2) if registros_originais > 0 else 0

                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("Registros Originais", registros_originais)
                        with col2:
                            st.metric("Registros Consolidados", registros_consolidados)
                        with col3:
                            st.metric("Redução", f"{reducao_percentual}%")

                        st.subheader("Preview dos Dados Consolidados")
                        st.dataframe(df_consolidado.head(10))

                        st.markdown("### Download")
                        nome_arquivo = f"{os.path.splitext(uploaded_file.name)[0]}_consolidado.csv"
                        st.markdown(get_csv_download_link(df_consolidado, nome_arquivo, "📥 Baixar Arquivo Consolidado"), unsafe_allow_html=True)

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
            """)

else:
    # Novo modo - Relatório de Variantes Genômicas
    with st.expander("ℹ️ Instruções e Informações", expanded=True):
        st.markdown("""
        Esta ferramenta gera um relatório compreensivo de variantes genômicas a partir de arquivos Excel.

        **Colunas Obrigatórias:**
        - Case_Record_ID: Identificador do caso/paciente
        - Gene: Nome do gene
        - CHR: Cromossomo
        - POS: Posição genômica
        - rsID: Identificador rs (dbSNP)
        - REF: Alelo de referência
        - ALT: Alelo alternativo
        - Genotype: Genótipo
        - Zygosity: Zigosidade (Homozygous/Heterozygous)
        - CLNSIG_Resolved: Significância clínica resolvida

        **Colunas Opcionais:**
        - CLNSIGCONF: Conflitos de interpretação
        - CLNDN: Condição clínica/doença

        **Como funciona:**
        1. Faça upload de um arquivo Excel (.xlsx ou .xls)
        2. O relatório será organizado por significância clínica (Pathogenic, Likely Pathogenic, VUS, etc.)
        3. Dentro de cada categoria, os dados serão agrupados por caso/paciente
        4. Para cada caso, as variantes serão organizadas por gene
        5. Variantes serão codificadas por cores baseadas na significância clínica
        """)

    uploaded_file = st.file_uploader("Selecione um arquivo Excel", type=["xlsx", "xls"])

    if uploaded_file is not None:
        try:
            # Ler arquivo Excel
            df = pd.read_excel(uploaded_file, engine='openpyxl')

            st.write(f"**Arquivo carregado:** {uploaded_file.name}")
            st.write(f"**Número de linhas:** {len(df)}")
            st.write(f"**Número de colunas:** {len(df.columns)}")

            # Processar dados genômicos
            df_genomic = parse_genomic_data(df)

            if df_genomic is not None:
                if st.button("🧬 Gerar Relatório de Variantes"):
                    with st.spinner("Gerando relatório compreensivo... Por favor, aguarde."):
                        generate_genomic_report(df_genomic)

        except Exception as e:
            st.error(f"Ocorreu um erro ao processar o arquivo: {str(e)}")
            st.info("""
            **Possíveis soluções:**
            - Verifique se o arquivo Excel contém todas as colunas obrigatórias
            - Certifique-se de que o arquivo está no formato Excel correto (.xlsx ou .xls)
            - Verifique se o arquivo não está corrompido
            - Verifique se os nomes das colunas estão corretos (Case_Record_ID, Gene, CHR, etc.)
            """)

# Rodapé
st.markdown("---")
st.markdown("Versão 2.0 - Desenvolvido para processamento e análise de dados genômicos")
