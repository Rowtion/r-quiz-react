import { useState } from 'react'
import { quizData } from './data/quizData'
import './App.css'

function App() {
  const [currentQuiz, setCurrentQuiz] = useState(null)
  const [currentQuestion, setCurrentQuestion] = useState(0)
  const [selectedOptions, setSelectedOptions] = useState({})
  const [showResults, setShowResults] = useState(false)
  const [nickname, setNickname] = useState('')
  const [isNicknameSet, setIsNicknameSet] = useState(false)

  const handleQuizSelect = (quizId) => {
    if (!isNicknameSet) return;
    setCurrentQuiz(quizId)
    setCurrentQuestion(0)
    setSelectedOptions({})
    setShowResults(false)
  }

  const handleNicknameSubmit = (e) => {
    e.preventDefault()
    if (nickname.trim()) {
      setIsNicknameSet(true)
    }
  }

  const handleOptionSelect = (questionIndex, optionIndex) => {
    if (!showResults) {
      setSelectedOptions({
        ...selectedOptions,
        [questionIndex]: optionIndex
      })
    }
  }

  const calculateScore = () => {
    let score = 0
    const quiz = quizData[currentQuiz]
    for (let i = 0; i < quiz.questions.length; i++) {
      // 只有当选项存在且等于正确答案时才计分
      if (selectedOptions[i] !== undefined && selectedOptions[i] === quiz.questions[i].correctAnswer) {
        score++
      }
    }
    return score
  }

  const getOptionClass = (questionIndex, optionIndex) => {
    if (!showResults) {
      return selectedOptions[questionIndex] === optionIndex ? 'option selected' : 'option'
    }

    const correctAnswer = quizData[currentQuiz].questions[questionIndex].correctAnswer
    if (optionIndex === correctAnswer) {
      return 'option correct'
    }
    if (selectedOptions[questionIndex] === optionIndex) {
      return 'option wrong'
    }
    if (selectedOptions[questionIndex] === undefined && optionIndex === correctAnswer) {
      return 'option correct unanswered'
    }
    return 'option'
  }

  const getOptionLabel = (index) => {
    return String.fromCharCode(65 + index); // 将0,1,2,3转换为A,B,C,D
  }

  const getExplanationClass = (questionIndex) => {
    if (!showResults) return 'explanation'
    const correctAnswer = quizData[currentQuiz].questions[questionIndex].correctAnswer
    if (selectedOptions[questionIndex] === undefined) {
      return 'explanation unanswered'
    }
    return selectedOptions[questionIndex] === correctAnswer 
      ? 'explanation correct'
      : 'explanation wrong'
  }

  const getQuestionStatusText = (questionIndex) => {
    if (!showResults) return null
    if (selectedOptions[questionIndex] === undefined) {
      return '(未作答)'
    }
    const selectedLabel = getOptionLabel(selectedOptions[questionIndex]);
    return selectedOptions[questionIndex] === quizData[currentQuiz].questions[questionIndex].correctAnswer
      ? `(正确 - 选择了${selectedLabel})`
      : `(错误 - 选择了${selectedLabel})`
  }

  return (
    <div className="app">
      {!isNicknameSet ? (
        <div className="nickname-container">
          <input
            type="text"
            placeholder="请输入你的昵称"
            value={nickname}
            onChange={(e) => setNickname(e.target.value)}
          />
          <button onClick={handleNicknameSubmit}>开始</button>
        </div>
      ) : !currentQuiz ? (
        <div className="quiz-selection">
          <h1>全代码搞定6+生信SCI分析与写作讲席营作业及测试</h1>
          <div className="student-info">
            欢迎你，<span className="nickname">{nickname}</span>
          </div>
          <div className="quiz-buttons">
            {Object.entries(quizData).map(([id, quiz]) => (
              <button key={id} onClick={() => handleQuizSelect(id)}>
                {quiz.title}
              </button>
            ))}
          </div>
        </div>
      ) : (
        <div className="quiz">
          {!showResults ? (
            <>
              <div className="quiz-header">
                <h2>{quizData[currentQuiz].title}</h2>
                <p>{quizData[currentQuiz].description}</p>
              </div>
              <div className="question">
                <div className="question-header">
                  {quizData[currentQuiz].questions[currentQuestion].category && (
                    <div className="category">
                      {quizData[currentQuiz].questions[currentQuestion].category}
                    </div>
                  )}
                  {quizData[currentQuiz].questions[currentQuestion].difficulty && (
                    <div className={`difficulty ${quizData[currentQuiz].questions[currentQuestion].difficulty}`}>
                      {quizData[currentQuiz].questions[currentQuestion].difficulty}
                    </div>
                  )}
                </div>
                <div>{quizData[currentQuiz].questions[currentQuestion].question}</div>
                {quizData[currentQuiz].questions[currentQuestion].code && (
                  <pre className="code-block">
                    <code>{quizData[currentQuiz].questions[currentQuestion].code}</code>
                  </pre>
                )}
                <div className="options">
                  {quizData[currentQuiz].questions[currentQuestion].options.map((option, optionIndex) => {
                    const isCorrect = showResults && optionIndex === quizData[currentQuiz].questions[currentQuestion].correctAnswer;
                    const isSelected = selectedOptions[currentQuestion] === optionIndex;
                    return (
                      <div
                        key={optionIndex}
                        className={`${getOptionClass(currentQuestion, optionIndex)} ${isCorrect ? 'correct-answer' : ''}`}
                        onClick={() => handleOptionSelect(currentQuestion, optionIndex)}
                      >
                        <span className="option-label">{getOptionLabel(optionIndex)}.</span>
                        {option}
                        {showResults && isCorrect && <span className="correct-indicator">✓</span>}
                      </div>
                    );
                  })}
                </div>
              </div>
              <div className="quiz-controls">
                <button
                  onClick={() => setCurrentQuestion(Math.max(0, currentQuestion - 1))}
                  disabled={currentQuestion === 0}
                >
                  上一题
                </button>
                {currentQuestion === quizData[currentQuiz].questions.length - 1 ? (
                  <button onClick={() => setShowResults(true)}>提交</button>
                ) : (
                  <button
                    onClick={() => setCurrentQuestion(Math.min(quizData[currentQuiz].questions.length - 1, currentQuestion + 1))}
                  >
                    下一题
                  </button>
                )}
              </div>
            </>
          ) : (
            <div className="results">
              <div className="quiz-result">
                <h2>测验结果</h2>
                <div className="student-info">
                  学生：<span className="nickname">{nickname}</span>
                </div>
                <div className="score-container">
                  <div className="score-label">你的得分</div>
                  <div className="score">{calculateScore()} / {quizData[currentQuiz].questions.length}</div>
                </div>
                <div className="score-breakdown">
                  <div className="score-item">
                    <div className="score-item-value correct-count">{calculateScore()}</div>
                    <div className="score-item-label">正确</div>
                  </div>
                  <div className="score-item">
                    <div className="score-item-value wrong-count">
                      {quizData[currentQuiz].questions.length - calculateScore()}
                    </div>
                    <div className="score-item-label">错误</div>
                  </div>
                </div>
                <div className="questions-overview">
                  {quizData[currentQuiz].questions.map((_, index) => {
                    const isCorrect = selectedOptions[index] === quizData[currentQuiz].questions[index].correctAnswer;
                    const isUnanswered = selectedOptions[index] === undefined;
                    return (
                      <div 
                        key={index} 
                        className={`question-dot ${isUnanswered ? 'unanswered' : (isCorrect ? 'correct' : 'wrong')}`}
                        title={`第 ${index + 1} 题: ${isUnanswered ? '未作答' : (isCorrect ? '正确' : '错误')}`}
                      >
                        {index + 1}
                      </div>
                    );
                  })}
                </div>
              </div>
            </div>
          )}
          {showResults && (
            <div className="return-button-container">
              <button onClick={() => {
                setCurrentQuiz(null);
                setShowResults(false);
                setSelectedOptions({});
              }}>
                返回选择测验
              </button>
              <div className="checkin-hint">
                请将测验结果截图上传至：
                <a href="https://w85j3acboi.feishu.cn/share/base/form/shrcnmgAOO4G1r9hAWn2TjDeoQd" 
                   target="_blank" 
                   rel="noopener noreferrer"
                   className="checkin-link">
                  飞书打卡链接
                </a>
              </div>
            </div>
          )}
        </div>
      )}
      {showResults && (
        <div className="questions-review">
          <h3>题目回顾</h3>
          {quizData[currentQuiz].questions.slice(0, -1).map((q, index) => (
            <div key={index} className="question-review">
              <div className="question-header">
                {q.category && <div className="category">{q.category}</div>}
                {q.difficulty && <div className={`difficulty ${q.difficulty}`}>{q.difficulty}</div>}
                <div className="question-status">{getQuestionStatusText(index)}</div>
              </div>
              <div>{q.question}</div>
              {q.code && <pre className="code-block">{q.code}</pre>}
              <div className="options">
                {q.options.map((option, optionIndex) => {
                  const isCorrect = showResults && optionIndex === q.correctAnswer;
                  return (
                    <div
                      key={optionIndex}
                      className={`${getOptionClass(index, optionIndex)} ${isCorrect ? 'correct-answer' : ''}`}
                    >
                      <span className="option-label">{getOptionLabel(optionIndex)}.</span>
                      {option}
                      {showResults && isCorrect && <span className="correct-indicator">✓</span>}
                    </div>
                  );
                })}
              </div>
              <div className={getExplanationClass(index)}>{q.explanation}</div>
            </div>
          ))}
        </div>
      )}
    </div>
  )
}

export default App
